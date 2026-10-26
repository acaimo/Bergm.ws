.draw_prior_phi <- function(prior, K, nsim) {
  U0 <- .chol_safe(prior$Sigma0)
  t(replicate(nsim, {
    mu <- prior$mu0 + .rmvn_u(U0)
    U <- .chol_safe(.riwish(prior$nu0, prior$S0))
    unlist(lapply(seq_len(K), function(k) mu + .rmvn_u(U)))
  }))
}

ws_predictive <- function(object, type = c("posterior", "prior"), prior = NULL,
                           nsim = 200, burnin = NULL, gof = TRUE, seed = NULL,
                           verbose = TRUE) {
  type <- match.arg(type)
  if (inherits(object, "bergm_w")) {
    model <- object$model
    if (is.null(prior)) prior <- object$prior
    aux <- object$settings$aux_iters
  } else if (inherits(object, "ws_model")) {
    if (type == "posterior")
      stop("the posterior predictive check needs a fit from bergm_w() or bergm_ws().", call. = FALSE)
    model <- object
    if (is.null(prior)) prior <- ws_prior(model)
    aux <- .default_aux_iters(model$data)
  } else stop("'object' must be a fit from bergm_w()/bergm_ws() or a ws_model.", call. = FALSE)
  if (is.null(burnin)) burnin <- 5L * aux
  if (!is.null(seed)) set.seed(seed)
  K <- model$K; p <- model$p; d <- K * p

  draws <- if (type == "posterior") {
    P <- do.call(rbind, lapply(object$chains, `[[`, "phi"))
    P[sample.int(nrow(P), nsim, replace = nrow(P) < nsim), , drop = FALSE]
  } else .draw_prior_phi(prior, K, nsim)

  sim <- .aux_sampler(model, as.integer(burnin))
  S <- matrix(NA_real_, nsim, d, dimnames = list(NULL, model$labels))
  G <- if (gof) matrix(NA_real_, nsim, length(.gof_stats(model$signed)) * K, dimnames = list(NULL, .gof_names(K, model$signed))) else NULL
  step <- max(1L, nsim %/% 5L)
  for (i in seq_len(nsim)) {
    s <- sim(draws[i, ], network = gof)
    S[i, ] <- s
    if (gof) G[i, ] <- .gof_simulated(attr(s, "network"), K, model$signed)
    if (verbose && i %% step == 0)
      message(sprintf("  %s predictive: %d / %d networks", type, i, nsim))
  }
  out <- list(type = type, model = model, nsim = nsim, burnin = burnin, draws = draws,
              stats = S, obs = model$s_obs_scaled,
              gof = G, gof_obs = if (gof) .gof_observed(model$data) else NULL)
  class(out) <- "ws_predictive"
  out
}

.mid_p <- function(S, obs)
  colMeans(sweep(S, 2, obs, "<")) + 0.5 * colMeans(sweep(S, 2, obs, "=="))

.pred_table <- function(S, obs) {
  q <- t(apply(S, 2, quantile, c(0.025, 0.975), names = FALSE))
  data.frame(observed = obs, mean = colMeans(S), `2.5%` = q[, 1], `97.5%` = q[, 2],
             `mid-p` = .mid_p(S, obs),
             outside = ifelse(obs < q[, 1] | obs > q[, 2], "*", ""),
             check.names = FALSE, row.names = colnames(S))
}

print.ws_predictive <- function(x, digits = 2, ...) {
  cat(if (x$type == "posterior") "Posterior" else "Prior",
      " predictive check: ", x$nsim, " simulated networks\n", sep = "")
  cat("(* = observed value outside the central 95% of the simulated values)\n\n")
  cat("Model statistics", if (any(x$model$scale != 1)) "(rescaled)" else "", "\n")
  t1 <- .pred_table(x$stats, x$obs)
  t1[1:4] <- lapply(t1[1:4], round, digits)
  t1[[5]] <- round(t1[[5]], 3)
  print(t1)
  if (!is.null(x$gof)) {
    cat("\nStatistics not in the model\n")
    t2 <- .pred_table(x$gof, x$gof_obs)
    t2[1:4] <- lapply(t2[1:4], round, digits)
    t2[[5]] <- round(t2[[5]], 3)
    print(t2)
  }
  invisible(x)
}

ws_aux_check <- function(fit, n = 100, factor = 20, seed = NULL) {
  if (!inherits(fit, "bergm_w")) stop("'fit' must come from bergm_w() or bergm_ws().", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  model <- fit$model
  aux <- fit$settings$aux_iters
  phi_bar <- colMeans(do.call(rbind, lapply(fit$chains, `[[`, "phi")))
  short <- .aux_sampler(model, aux)
  long <- .aux_sampler(model, as.integer(aux * factor))
  Ss <- t(replicate(n, short(phi_bar)))
  Sl <- t(replicate(n, long(phi_bar)))
  sd_l <- apply(Sl, 2, sd)
  sd_l[sd_l == 0] <- 1
  diff <- (colMeans(Ss) - colMeans(Sl)) / sd_l
  se <- sqrt(apply(Ss, 2, var) / n + apply(Sl, 2, var) / n) / sd_l
  out <- data.frame(statistic = model$labels, mean_short = colMeans(Ss),
                    mean_long = colMeans(Sl), sd_long = sd_l, std_diff = diff,
                    z = diff / pmax(se, 1e-8), row.names = NULL)
  attr(out, "aux_iters") <- aux
  attr(out, "factor") <- factor
  class(out) <- c("ws_aux_check", "data.frame")
  out
}

print.ws_aux_check <- function(x, ...) {
  cat("Auxiliary chains of ", attr(x, "aux_iters"), " steps vs ", attr(x, "factor"),
      " times longer, at the posterior mean\n", sep = "")
  zc <- stats::qnorm(1 - 0.025 / nrow(x))
  cat("(std_diff = difference of means in units of the long-chain sd; flagged if\n",
      " |std_diff| > 0.25 and |z| > ", sprintf("%.2f", zc), ", Bonferroni over ", nrow(x),
      " statistics)\n\n", sep = "")
  y <- as.data.frame(unclass(x))
  y[-1] <- lapply(y[-1], round, 3)
  print(y, row.names = FALSE)
  flag <- abs(x$std_diff) > 0.25 & abs(x$z) > zc
  if (any(flag)) {
    cat("\nShort chains differ from long ones for:", paste(x$statistic[flag], collapse = ", "),
        "\nConsider increasing 'aux_iters'.\n")
  } else cat("\nNo appreciable difference: the auxiliary chains are long enough.\n")
  invisible(x)
}
