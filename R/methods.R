as.mcmc.list.bergm_w <- function(x, what = c("phi", "mu", "sd"), ...) {
  what <- match.arg(what)
  m <- x$model
  st <- x$settings
  lst <- lapply(x$chains, function(ch) {
    M <- switch(what,
      phi = ch$phi,
      mu = ch$mu,
      sd = do.call(rbind, lapply(seq_len(nrow(ch$sigma)), function(i)
        sqrt(diag(.unvech(ch$sigma[i, ], m$p))))))
    M <- matrix(M, ncol = if (what == "phi") m$K * m$p else m$p)
    colnames(M) <- switch(what, phi = m$labels, mu = paste0("mu:", m$terms),
                          sd = paste0("sd:", m$terms))
    mcmc(M, start = st$burnin + st$thin, thin = st$thin)
  })
  mcmc.list(lst)
}

.post_table <- function(ml, probs) {
  M <- do.call(rbind, lapply(ml, as.matrix))
  q <- t(apply(M, 2, quantile, probs = probs, names = FALSE))
  colnames(q) <- paste0(format(100 * probs, trim = TRUE), "%")
  res <- data.frame(mean = colMeans(M), sd = apply(M, 2, sd), q, check.names = FALSE)
  res$Rhat <- if (length(ml) > 1) tryCatch(
    gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1],
    error = function(e) rep(NA_real_, ncol(M))) else NA_real_
  res$ESS <- effectiveSize(ml)
  res
}

summary.bergm_w <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  m <- object$model
  K <- m$K; p <- m$p
  phi <- as.mcmc.list(object, "phi")
  burn <- seq_len(object$settings$burnin)
  acc <- vapply(object$chains, function(ch) mean(ch$accept[-burn, ]), numeric(1))
  P <- do.call(rbind, lapply(phi, as.matrix))
  C <- suppressWarnings(stats::cor(P))
  cross <- outer(rep(seq_len(K), each = p), rep(seq_len(K), each = p), "!=")
  sim_share <- mean(vapply(object$chains, function(ch)
    if (is.null(ch$n_aux)) 1 else mean(ch$n_aux[-burn, ]), numeric(1)))
  out <- list(
    sim_share = sim_share, delayed = isTRUE(object$settings$delayed),
    mu = .post_table(as.mcmc.list(object, "mu"), probs),
    phi = .post_table(phi, probs),
    sd = .post_table(as.mcmc.list(object, "sd"), probs),
    acceptance = acc,
    max_cross_corr = if (K > 1) max(abs(C[cross]), na.rm = TRUE) else NA_real_,
    corr = C, model = m, update = object$update, n_chains = object$n_chains,
    draws = nrow(P), settings = object$settings, elapsed = object$elapsed)
  class(out) <- "summary.bergm_w"
  out
}

.print_post <- function(tab, digits = 3) {
  num <- vapply(tab, is.numeric, logical(1))
  f <- tab
  f[num] <- lapply(names(tab)[num], function(nm) {
    v <- tab[[nm]]
    if (nm == "ESS") formatC(round(v), format = "d", big.mark = ",") else
      if (nm == "Rhat") .fmt(v, 3) else .fmt(v, digits)
  })
  print(f, right = TRUE)
}

print.summary.bergm_w <- function(x, digits = 3, ...) {
  m <- x$model
  cat("Hierarchical multi-layer dissolution ERGM, ", if (m$signed) "signed" else "unsigned",
      " weighted network\n  adaptive exchange (", x$update, " update)\n", sep = "")
  cat(sprintf("  %d chains, %s retained draws, %d auxiliary MCMC steps per exchange, %.0f s\n",
              x$n_chains, format(x$draws, big.mark = ","), x$settings$aux_iters, x$elapsed))
  cat("  acceptance after burn-in:", paste(.fmt(x$acceptance, 2), collapse = ", "), "\n")
  if (x$delayed)
    cat(sprintf("  delayed acceptance: %.0f%% of the proposals needed an auxiliary network\n",
                100 * x$sim_share))
  if (!is.null(.scale_note(m))) cat(" ", .scale_note(m), "\n")
  cat("\nHierarchical means\n")
  .print_post(x$mu, digits)
  for (k in seq_len(m$K)) {
    cat("\nLayer ", k, " (x_", k - 1, " -> x_", k, ")\n", sep = "")
    .print_post(x$phi[(k - 1) * m$p + seq_len(m$p), ], digits)
  }
  cat("\nBetween-layer standard deviations, sqrt(diag(Sigma))\n")
  .print_post(x$sd, digits)
  r <- c(x$phi$Rhat, x$mu$Rhat)
  e <- c(x$phi$ESS, x$mu$ESS)
  cat("\nConvergence: max Rhat ", .fmt(max(r, na.rm = TRUE), 3),
      ", min ESS ", round(min(e)),
      if (!is.na(x$max_cross_corr)) paste0(", max |corr| between layers ", .fmt(x$max_cross_corr, 2)),
      "\n", sep = "")
  if (any(r > 1.05, na.rm = TRUE) || any(e < 100))
    cat("  Some parameters have Rhat > 1.05 or ESS < 100: run longer chains.\n")
  invisible(x)
}

print.bergm_w <- function(x, ...) {
  s <- summary(x)
  cat(if (x$model$signed) "bergm_ws" else "bergm_w", ": ", x$n_chains, " chains (", x$update, " update), ",
      format(s$draws, big.mark = ","), " retained draws\n", sep = "")
  cat("Posterior means of the layer-specific parameters:\n")
  print(round(coef(x), 3))
  cat("max Rhat", .fmt(max(c(s$phi$Rhat, s$mu$Rhat), na.rm = TRUE), 3),
      "| min ESS", round(min(c(s$phi$ESS, s$mu$ESS))), "\n")
  cat("Use summary() for intervals and diagnostics, as.mcmc.list() for the draws.\n")
  invisible(x)
}

coef.bergm_w <- function(object, what = c("phi", "mu"), scale = c("model", "ergm"), ...) {
  what <- match.arg(what)
  scale <- match.arg(scale)
  m <- object$model
  if (what == "mu") {
    v <- colMeans(do.call(rbind, lapply(object$chains, `[[`, "mu")))
    names(v) <- m$terms
    return(v)
  }
  v <- colMeans(do.call(rbind, lapply(object$chains, `[[`, "phi")))
  if (scale == "ergm") v <- v * m$scale
  matrix(v, m$K, m$p, byrow = TRUE,
         dimnames = list(paste0("layer ", seq_len(m$K)), m$terms))
}
