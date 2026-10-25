.default_aux_iters <- function(data) {
  E <- vapply(data$x, network::network.edgecount, numeric(1))
  n_vary <- (if (data$signed) 2 * E[1] else E[1]) + sum(E[seq_len(data$K - 1) + 1])
  as.integer(max(2000, 10 * n_vary))
}

.aux_sampler <- function(model, aux_iters) {
  ctrl <- control.simulate.formula(MCMC.burnin = aux_iters)
  state <- simulate(model$formula, coef = NULL, basis = model$data$y,
                    constraints = model$data$constraint, nsim = 1,
                    output = "stats", control = ctrl,
                    return.args = "ergm_state")$object
  ctrl$MCMC.samplesize <- 1
  sc <- unname(model$scale)
  function(phi, network = FALSE) {
    phi <- unname(phi)
    res <- ergm_MCMC_sample(state, ctrl, theta = phi * sc)
    s <- .stats_vec(res) * sc
    if (network) attr(s, "network") <- as.network(res$networks[[1]])
    s
  }
}

.lhs_formula <- function(model) {
  env <- new.env(parent = environment(model$formula))
  assign("y", model$data$y, envir = env)
  as.formula(paste("y ~", model$rhs), env = env)
}

.mple_design <- function(model) {
  mp <- tryCatch(suppressMessages(suppressWarnings(
    ergmMPLE(.lhs_formula(model), constraints = model$data$constraint, output = "matrix"))),
    error = function(e) NULL)
  if (is.null(mp) || ncol(mp$predictor) != model$K * model$p) return(NULL)
  mp
}

.mple_map <- function(model, prior, maxit = 100) {
  mp <- .mple_design(model)
  if (is.null(mp)) return(NULL)
  K <- model$K
  X <- sweep(mp$predictor, 2, unname(model$scale), "*")
  y <- mp$response
  w <- mp$weights
  P <- chol2inv(.chol_safe(.prior_cov_phi(prior, K)))
  m0 <- rep(prior$mu0, K)
  obj <- function(phi) {
    eta <- drop(X %*% phi)
    sum(w * (y * eta - pmax(eta, 0) - log1p(exp(-abs(eta))))) -
      0.5 * sum((phi - m0) * (P %*% (phi - m0)))
  }
  phi <- m0
  f0 <- obj(phi)
  for (it in seq_len(maxit)) {
    pr <- plogis(drop(X %*% phi))
    g <- drop(crossprod(X, w * (y - pr))) - drop(P %*% (phi - m0))
    H <- crossprod(X * (w * pr * (1 - pr)), X) + P
    delta <- drop(chol2inv(.chol_safe(H)) %*% g)
    step <- 1
    repeat {
      f1 <- obj(phi + step * delta)
      if (f1 >= f0 - 1e-10 || step < 1e-4) break
      step <- step / 2
    }
    phi <- phi + step * delta
    if (abs(f1 - f0) < 1e-10 * (1 + abs(f0))) break
    f0 <- f1
  }
  phi
}

.sim_stats <- function(model, phi, nsim, aux_iters) {
  sims <- tryCatch(suppressMessages(as.matrix(simulate(
    model$formula, coef = unname(phi) * unname(model$scale), basis = model$data$y,
    constraints = model$data$constraint, nsim = nsim, output = "stats",
    control = control.simulate.formula(MCMC.burnin = 10L * aux_iters,
                                       MCMC.interval = max(1L, aux_iters %/% 4L))))),
    error = function(e) NULL)
  if (is.null(sims)) return(NULL)
  sweep(sims, 2, unname(model$scale), "*")
}

.ws_init <- function(model, prior, init, prop_cov, aux_iters, n_info, n_newton) {
  K <- model$K; p <- model$p; d <- K * p
  s_obs <- unname(model$s_obs_scaled)
  phi <- rep(0, d)
  start_type <- "zero"
  if (is.numeric(init)) {
    phi <- rep_len(as.numeric(init), d)
    start_type <- "user"
  } else if (identical(init, "mple")) {
    pm <- .mple_map(model, prior)
    if (!is.null(pm) && all(is.finite(pm))) {
      phi <- pm
      start_type <- "penalised MPLE"
    } else {
      warning("the pseudo-likelihood could not be computed; chains start at zero.",
              call. = FALSE)
    }
  } else if (!identical(init, "zero")) {
    stop("'init' must be \"mple\", \"zero\" or a numeric vector.", call. = FALSE)
  }

  P <- chol2inv(.chol_safe(.prior_cov_phi(prior, K)))
  m0 <- rep(prior$mu0, K)
  B <- NULL
  S <- NULL
  newton_done <- 0L
  if (identical(prop_cov, "information") || n_newton > 0) {
    for (step in 0:n_newton) {
      S <- .sim_stats(model, phi, n_info, aux_iters)
      if (is.null(S) || any(apply(S, 2, sd) == 0)) break
      I <- cov(S)
      if (step == n_newton) break
      H <- I + P
      g <- (s_obs - colMeans(S)) - drop(P %*% (phi - m0))
      delta <- drop(chol2inv(.chol_safe(H)) %*% g)
      len <- sqrt(sum(delta * drop(H %*% delta)))
      if (len > 3 * sqrt(d)) delta <- delta * 3 * sqrt(d) / len
      phi <- phi + delta
      newton_done <- newton_done + 1L
    }
  }
  surrogate <- NULL
  if (!is.null(S) && all(apply(S, 2, sd) > 0))
    surrogate <- list(center = phi, g = s_obs - colMeans(S), H = cov(S))
  fit_check <- NULL
  if (!is.null(S)) {
    sds <- apply(S, 2, sd)
    fit_check <- data.frame(statistic = model$labels, observed = s_obs,
                            simulated_mean = colMeans(S), simulated_sd = sds,
                            z = (s_obs - colMeans(S)) / pmax(sds, 1e-8), row.names = NULL)
  }
  if (is.matrix(prop_cov)) {
    if (any(dim(prop_cov) != d)) stop("'prop_cov' must be ", d, " x ", d, ".", call. = FALSE)
    B <- prop_cov
  } else if (identical(prop_cov, "information")) {
    if (!is.null(S) && all(apply(S, 2, sd) > 0)) {
      B <- chol2inv(.chol_safe(cov(S) + P))
    } else {
      warning("some statistics do not vary in networks simulated at the starting values; ",
              "using a diagonal initial proposal.", call. = FALSE)
    }
  } else if (!identical(prop_cov, "diagonal")) {
    stop("'prop_cov' must be \"information\", \"diagonal\" or a matrix.", call. = FALSE)
  }
  from_info <- !is.null(B)
  if (is.null(B)) B <- diag(0.01, d)
  dimnames(B) <- list(model$labels, model$labels)
  list(phi = setNames(phi, model$labels), B = B, start = start_type,
       newton_steps = newton_done, from_info = from_info, fit_check = fit_check,
       surrogate = surrogate)
}

.blocks <- function(update, K, p) {
  if (update == "joint") list(seq_len(K * p)) else
    lapply(seq_len(K), function(k) (k - 1L) * p + seq_len(p))
}

.init_state <- function(phi, prior, K, p, update, B0, from_info) {
  blocks <- .blocks(update, K, p)
  list(phi = unname(phi), mu = colMeans(matrix(phi, K, p, byrow = TRUE)),
       Sigma = prior$E_Sigma, gam = 2.38 / sqrt(2 * lengths(blocks)),
       shape = lapply(blocks, function(ix) B0[ix, ix, drop = FALSE]),
       posterior_shape = rep(from_info, length(blocks)))
}

.ws_run <- function(model, prior, update, state, B0, from_info, n_iter, store_from,
                     thin, aux_iters, adapt_cov, adapt_start, adapt_stop, adapt_every,
                     adapt_scaling, target_accept, seed, surrogate = NULL,
                     collapse_mu = TRUE, verbose = FALSE, id = 1L, it_offset = 0L,
                     it_total = n_iter) {
  set.seed(seed)
  K <- model$K; p <- model$p; d <- K * p
  s_obs <- unname(model$s_obs_scaled)
  aux <- .aux_sampler(model, aux_iters)

  blocks <- .blocks(update, K, p)
  nb <- length(blocks)
  bl <- lengths(blocks)
  opt <- 2.38 / sqrt(2 * bl)
  shape0 <- lapply(blocks, function(ix) B0[ix, ix, drop = FALSE])
  idx_k <- lapply(seq_len(K), function(k) (k - 1L) * p + seq_len(p))
  blk_layers <- lapply(blocks, function(ix) unique((ix - 1L) %/% p + 1L))

  phi <- state$phi; mu <- state$mu; Sigma <- state$Sigma
  gam <- state$gam; shape <- state$shape; posterior_shape <- state$posterior_shape
  U_S <- .chol_safe(Sigma)
  U_prop <- lapply(seq_len(nb), function(b) .chol_safe(gam[b]^2 * shape[[b]]))
  m_marg <- rep(prior$mu0, K)
  VJ <- kronecker(matrix(1, K, K), prior$Sigma0)
  marg_chol <- function(S) .chol_safe(kronecker(diag(K), S) + VJ)
  U_V <- if (collapse_mu) marg_chol(Sigma) else NULL

  n_save <- max(0L, (n_iter - store_from) %/% thin)
  out_phi <- matrix(NA_real_, n_save, d)
  out_mu <- matrix(NA_real_, n_save, p)
  out_sig <- matrix(NA_real_, n_save, p * (p + 1) / 2)
  acc <- matrix(0L, n_iter, nb)
  n_aux <- matrix(0L, n_iter, nb)
  n_hist <- if (adapt_cov) min(adapt_stop, n_iter) else 0L
  hist <- if (n_hist > 0) matrix(NA_real_, n_hist, d) else NULL

  delayed <- !is.null(surrogate)
  sur <- if (delayed) {
    cc <- surrogate$center; gg <- surrogate$g; HH <- surrogate$H
    function(x) { u <- x - cc; sum(gg * u) - 0.5 * sum(u * drop(HH %*% u)) }
  } else NULL
  sur_phi <- if (delayed) sur(phi) else 0
  report <- if (verbose) unique(round(seq(it_total / 10, it_total, length.out = 10))) else integer(0)

  for (it in seq_len(n_iter)) {
    for (b in seq_len(nb)) {
      ix <- blocks[[b]]
      prop <- phi
      prop[ix] <- phi[ix] + .rmvn_u(U_prop[[b]])
      if (collapse_mu) {
        lprior <- .ldmvn_u(prop, m_marg, U_V) - .ldmvn_u(phi, m_marg, U_V)
      } else {
        lprior <- 0
        for (k in blk_layers[[b]]) {
          kk <- idx_k[[k]]
          lprior <- lprior + .ldmvn_u(prop[kk], mu, U_S) - .ldmvn_u(phi[kk], mu, U_S)
        }
      }
      if (delayed) {
        sur_prop <- sur(prop)
        dsur <- sur_prop - sur_phi
        if (!(log(runif(1)) < dsur + lprior)) next
        s_aux <- aux(prop)
        n_aux[it, b] <- 1L
        la <- sum((phi - prop) * (s_aux - s_obs)) - dsur
      } else {
        s_aux <- aux(prop)
        n_aux[it, b] <- 1L
        la <- sum((phi - prop) * (s_aux - s_obs)) + lprior
      }
      if (is.finite(la) && log(runif(1)) < la) {
        phi <- prop
        if (delayed) sur_phi <- sur_prop
        acc[it, b] <- 1L
      }
    }
    if (!is.null(hist) && it <= n_hist) hist[it, ] <- phi

    if (it >= adapt_start && it <= adapt_stop && it %% adapt_every == 0) {
      rate <- colMeans(acc[(it - adapt_every + 1L):it, , drop = FALSE])
      for (b in seq_len(nb)) {
        if (!is.null(hist)) {
          h <- hist[ceiling(it / 2):it, blocks[[b]], drop = FALSE]
          if (nrow(h) >= 10L * bl[b]) {
            C <- cov(h)
            if (all(is.finite(C)) && sum(diag(C)) > 0) {
              base <- if (from_info) shape0[[b]] else
                shape0[[b]] * sum(diag(C)) / sum(diag(shape0[[b]]))
              shape[[b]] <- 0.95 * C + 0.05 * base
              if (!posterior_shape[b]) {
                gam[b] <- opt[b]
                posterior_shape[b] <- TRUE
              }
            }
          }
        }
        gam[b] <- gam[b] * (if (rate[b] > target_accept) 1 + adapt_scaling else 1 - adapt_scaling)
        gam[b] <- min(max(gam[b], 1e-3), 10)
        U_prop[[b]] <- .chol_safe(gam[b]^2 * shape[[b]])
      }
    }

    Phi <- matrix(phi, K, p, byrow = TRUE)
    Sinv <- chol2inv(U_S)
    Smu <- chol2inv(.chol_safe(prior$Sigma0_inv + K * Sinv))
    m <- Smu %*% (prior$Sigma0_inv %*% prior$mu0 + K * Sinv %*% colMeans(Phi))
    mu <- drop(m) + .rmvn_u(.chol_safe(Smu))
    Sigma <- .riwish(prior$nu0 + K, prior$S0 + crossprod(sweep(Phi, 2, mu)))
    U_S <- .chol_safe(Sigma)
    if (collapse_mu) U_V <- marg_chol(Sigma)

    if (it > store_from && (it - store_from) %% thin == 0) {
      j <- (it - store_from) %/% thin
      out_phi[j, ] <- phi
      out_mu[j, ] <- mu
      out_sig[j, ] <- .vech(Sigma)
    }
    if ((it + it_offset) %in% report)
      message(sprintf("  chain %d: %3.0f%%", id, 100 * (it + it_offset) / it_total))
  }
  list(phi = out_phi, mu = out_mu, sigma = out_sig, accept = acc, n_aux = n_aux,
       state = list(phi = phi, mu = mu, Sigma = Sigma, gam = gam, shape = shape,
                    posterior_shape = posterior_shape),
       hist = if (!is.null(hist)) hist[ceiling(n_hist / 2):n_hist, , drop = FALSE] else NULL)
}

.stop_if_failed <- function(res) {
  bad <- vapply(res, function(r) inherits(r, "error") || inherits(r, "try-error") ||
                  is.null(r$state), logical(1))
  if (any(bad)) {
    msg <- res[[which(bad)[1]]]
    stop("chain ", which(bad)[1], " failed: ",
         if (inherits(msg, "error")) conditionMessage(msg) else as.character(msg), call. = FALSE)
  }
  invisible(TRUE)
}

.bergm_fit <- function(model, prior = ws_prior(model), update = c("joint", "layer"),
                    n_iter = 10000, burnin = floor(n_iter / 4), thin = 1,
                    n_chains = 4, cores = 1, aux_iters = NULL, init = "mple",
                    prop_cov = "information", adapt_cov = TRUE, pool_adaptation = TRUE,
                    adapt_every = NULL, adapt_start = NULL, adapt_scaling = 0.1,
                    target_accept = 0.234, delayed = TRUE, collapse_mu = TRUE,
                    n_info = 200, n_newton = 2,
                    dispersion = 4, seed = NULL, verbose = TRUE) {
  if (!inherits(model, "ws_model")) stop("'model' must be created with ws_model().", call. = FALSE)
  if (!inherits(prior, "ws_prior")) stop("'prior' must be created with ws_prior().", call. = FALSE)
  if (prior$p != model$p) stop("the prior has ", prior$p, " statistics per layer, the model ",
                               model$p, ".", call. = FALSE)
  update <- match.arg(update)
  K <- model$K; p <- model$p
  n_iter <- as.integer(n_iter); burnin <- as.integer(burnin); thin <- as.integer(thin)
  if (burnin >= n_iter) stop("'burnin' must be smaller than 'n_iter'.", call. = FALSE)
  if ((n_iter - burnin) %/% thin < 2) stop("too few retained draws.", call. = FALSE)
  n_chains <- as.integer(n_chains)
  if (is.null(aux_iters)) aux_iters <- .default_aux_iters(model$data)
  aux_iters <- as.integer(aux_iters)
  if (is.null(adapt_every)) adapt_every <- max(20L, burnin %/% 20L)
  if (is.null(adapt_start)) adapt_start <- 2L * adapt_every
  adapt_every <- as.integer(adapt_every); adapt_start <- as.integer(adapt_start)
  if (burnin < adapt_start + adapt_every)
    warning("the burn-in is too short for the adaptation to take place; ",
            "increase 'burnin' or reduce 'adapt_every'.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  t0 <- proc.time()[["elapsed"]]

  if (verbose) message("Starting values and initial proposal covariance ...")
  ini <- .ws_init(model, prior, init, prop_cov, aux_iters, n_info, n_newton)
  if (delayed && is.null(ini$surrogate)) {
    warning("no Fisher information at the starting value: delayed acceptance switched off.",
            call. = FALSE)
    delayed <- FALSE
  }
  sur <- if (delayed) ini$surrogate else NULL
  U0 <- .chol_safe(ini$B)
  starts <- lapply(seq_len(n_chains), function(i)
    if (i == 1L) ini$phi else ini$phi + sqrt(dispersion) * .rmvn_u(U0))
  states <- lapply(starts, .init_state, prior = prior, K = K, p = p, update = update,
                   B0 = ini$B, from_info = ini$from_info)
  chain_seeds <- sample.int(.Machine$integer.max - 1L, n_chains)

  par_ok <- cores > 1 && n_chains > 1 && .Platform$OS.type != "windows"
  run_all <- function(FUN) {
    res <- if (par_ok) {
      mclapply(seq_len(n_chains), function(i) tryCatch(FUN(i), error = function(e) e),
               mc.cores = min(cores, n_chains), mc.set.seed = FALSE)
    } else lapply(seq_len(n_chains), FUN)
    .stop_if_failed(res)
    res
  }
  pool <- isTRUE(pool_adaptation) && adapt_cov && n_chains > 1 &&
    floor(0.75 * burnin) >= adapt_start + adapt_every
  if (verbose)
    message(sprintf("Running %d chain(s) x %d iterations (%s update%s, %d auxiliary MCMC steps)%s ...",
                    n_chains, n_iter, update, if (delayed) ", delayed acceptance" else "",
                    aux_iters, if (par_ok) sprintf(" on %d cores", min(cores, n_chains)) else ""))
  vb <- verbose && !par_ok

  run_phase <- function(st, n, store_from, adapt_cov, adapt_start, adapt_stop,
                        from_info, seed_offset, it_offset)
    run_all(function(i)
      .ws_run(model, prior, update, st[[i]], ini$B, from_info, n_iter = n,
               store_from = store_from, thin = thin, aux_iters = aux_iters,
               adapt_cov = adapt_cov, adapt_start = adapt_start, adapt_stop = adapt_stop,
               adapt_every = adapt_every, adapt_scaling = adapt_scaling,
               target_accept = target_accept, seed = chain_seeds[i] + seed_offset,
               surrogate = sur, collapse_mu = collapse_mu, verbose = vb, id = i,
               it_offset = it_offset, it_total = n_iter))

  if (pool) {
    b1 <- as.integer(floor(0.75 * burnin))
    A <- run_phase(states, b1, b1, TRUE, adapt_start, b1, ini$from_info, 0L, 0L)
    blocks <- .blocks(update, K, p)
    pooled <- lapply(blocks, function(ix) {
      C <- Reduce(`+`, lapply(A, function(r) cov(r$hist[, ix, drop = FALSE]))) / length(A)
      0.95 * C + 0.05 * ini$B[ix, ix, drop = FALSE]
    })
    gam_A <- exp(colMeans(do.call(rbind, lapply(A, function(r) log(r$state$gam)))))
    states <- lapply(A, function(r) {
      st <- r$state
      st$shape <- pooled
      st$gam <- gam_A
      st$posterior_shape[] <- TRUE
      st
    })
    b2 <- burnin - b1
    B <- run_phase(states, b2, b2, FALSE, adapt_every, b2, TRUE, 1L, b1)
    gam_B <- exp(colMeans(do.call(rbind, lapply(B, function(r) log(r$state$gam)))))
    states <- lapply(B, function(r) {
      st <- r$state
      st$gam <- gam_B
      st
    })
    C3 <- run_phase(states, n_iter - burnin, 0L, FALSE, 1L, 0L, TRUE, 2L, burnin)
    runs <- lapply(seq_len(n_chains), function(i)
      list(C3[[i]], accept = rbind(A[[i]]$accept, B[[i]]$accept, C3[[i]]$accept),
           n_aux = rbind(A[[i]]$n_aux, B[[i]]$n_aux, C3[[i]]$n_aux)))
  } else {
    R <- run_phase(states, n_iter, burnin, adapt_cov, adapt_start, burnin,
                   ini$from_info, 0L, 0L)
    runs <- lapply(R, function(r) list(r, accept = r$accept, n_aux = r$n_aux))
  }
  chains <- lapply(seq_len(n_chains), function(i) {
    r <- runs[[i]][[1]]
    list(phi = r$phi, mu = r$mu, sigma = r$sigma, accept = runs[[i]]$accept,
         n_aux = runs[[i]]$n_aux, gamma = r$state$gam, shape = r$state$shape,
         start = unname(starts[[i]]), seed = chain_seeds[i])
  })

  out <- list(model = model, prior = prior, update = update, chains = chains,
              n_chains = n_chains,
              settings = list(n_iter = n_iter, burnin = burnin, thin = thin,
                              aux_iters = aux_iters, adapt_cov = adapt_cov,
                              pool_adaptation = pool,
                              adapt_start = adapt_start, adapt_every = adapt_every,
                              adapt_scaling = adapt_scaling, target_accept = target_accept,
                              delayed = delayed, collapse_mu = collapse_mu,
                              n_info = n_info, n_newton = n_newton,
                              dispersion = dispersion, cores = cores, seed = seed),
              init = ini, elapsed = proc.time()[["elapsed"]] - t0, call = match.call())
  class(out) <- c(if (model$signed) "bergm_ws", "bergm_w")
  if (verbose) message(sprintf("Done in %.1f s.", out$elapsed))
  out
}

.as_model <- function(data, formula, scale, labels, signed) {
  if (inherits(data, "ws_model")) {
    if (signed && !data$signed) stop("the model is for an unsigned network: use bergm_w().", call. = FALSE)
    if (!signed && data$signed) stop("the model is for a signed network: use bergm_ws().", call. = FALSE)
    if (!is.null(formula)) warning("'data' is a ws_model: 'formula' is ignored.", call. = FALSE)
    return(data)
  }
  if (!inherits(data, "ws_data"))
    stop("'data' must be created with ws_data() or ws_model().", call. = FALSE)
  if (signed && !data$signed)
    stop("the network has no signs: use bergm_w(), or give 'signs' to ws_data().", call. = FALSE)
  if (!signed && data$signed) data <- ws_data(data$weights, data$thresholds, attributes = data$attributes)
  if (is.null(formula)) stop("'formula' is required.", call. = FALSE)
  ws_model(data, formula = formula, scale = scale, labels = labels)
}

.run_bergm <- function(model, prior, env, call) {
  a <- mget(setdiff(names(formals(.bergm_fit)), c("model", "prior")), envir = env)
  fit <- do.call(.bergm_fit, c(list(model, if (is.null(prior)) ws_prior(model) else prior), a))
  fit$call <- call
  fit
}

bergm_w <- function(data, formula = NULL, scale = NULL, labels = NULL, prior = NULL,
                    update = c("joint", "layer"), n_iter = 10000,
                    burnin = floor(n_iter / 4), thin = 1, n_chains = 4, cores = 1,
                    aux_iters = NULL, init = "mple", prop_cov = "information",
                    adapt_cov = TRUE, pool_adaptation = TRUE, adapt_every = NULL,
                    adapt_start = NULL, adapt_scaling = 0.1, target_accept = 0.234,
                    delayed = TRUE, collapse_mu = TRUE, n_info = 200, n_newton = 2,
                    dispersion = 4, seed = NULL, verbose = TRUE) {
  model <- .as_model(data, formula, scale, labels, signed = FALSE)
  .run_bergm(model, prior, environment(), match.call())
}

bergm_ws <- function(data, formula = NULL, scale = NULL, labels = NULL, prior = NULL,
                     update = c("joint", "layer"), n_iter = 10000,
                     burnin = floor(n_iter / 4), thin = 1, n_chains = 4, cores = 1,
                     aux_iters = NULL, init = "mple", prop_cov = "information",
                     adapt_cov = TRUE, pool_adaptation = TRUE, adapt_every = NULL,
                     adapt_start = NULL, adapt_scaling = 0.1, target_accept = 0.234,
                     delayed = TRUE, collapse_mu = TRUE, n_info = 200, n_newton = 2,
                     dispersion = 4, seed = NULL, verbose = TRUE) {
  model <- .as_model(data, formula, scale, labels, signed = TRUE)
  .run_bergm(model, prior, environment(), match.call())
}
