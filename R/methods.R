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
  cat("Use summary() for intervals and diagnostics, plot() for traces.\n")
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

.plot_intervals <- function(x, main = NULL) {
  pal <- .ws_pal
  m <- x$model
  K <- m$K; p <- m$p
  phi <- do.call(rbind, lapply(x$chains, `[[`, "phi"))
  mu <- do.call(rbind, lapply(x$chains, `[[`, "mu"))
  qs <- function(v) quantile(v, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
  lay_col <- colorRampPalette(c("#6da7ec", "#2a78d6", "#0d366b"), space = "Lab")(max(K, 2))
  if (K == 1) lay_col <- "#2a78d6"
  cols <- c(lay_col[seq_len(K)], pal$ink)
  nsub <- K + 1
  rows <- list()
  for (j in seq_len(p)) {
    for (k in seq_len(K)) rows[[length(rows) + 1]] <- list(j = j, s = k, q = qs(phi[, (k - 1) * p + j]))
    rows[[length(rows) + 1]] <- list(j = j, s = nsub, q = qs(mu[, j]))
  }
  y <- vapply(rows, function(r) -((r$j - 1) * (nsub + 1.2) + r$s), numeric(1))
  xr <- range(unlist(lapply(rows, `[[`, "q")), 0)
  xr <- xr + c(-1, 1) * 0.04 * diff(xr)
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  .ws_par(mar = c(5.2, 10, 3.6, 1.2))
  plot.new()
  plot.window(xlim = xr, ylim = c(min(y) - 0.6, max(y) + 0.6))
  abline(v = pretty(xr), col = pal$grid, lwd = 0.7)
  abline(v = 0, col = pal$axis, lwd = 1.1)
  for (i in seq_along(rows)) {
    r <- rows[[i]]
    segments(r$q[1], y[i], r$q[5], y[i], col = cols[r$s], lwd = 1.1, lend = 1)
    segments(r$q[2], y[i], r$q[4], y[i], col = cols[r$s], lwd = 4.2, lend = 1)
    points(r$q[3], y[i], pch = if (r$s == nsub) 23 else 21, bg = cols[r$s],
           col = pal$surface, cex = if (r$s == nsub) 1.25 else 1.05, lwd = 1.3)
  }
  axis(1, at = pretty(xr), lwd = 0, lwd.ticks = 0.6, cex.axis = 0.75)
  centre <- vapply(seq_len(p), function(j) mean(y[vapply(rows, `[[`, 0, "j") == j]), numeric(1))
  mtext(m$terms, side = 2, at = centre, line = 0.6, las = 1, cex = 0.8, col = pal$ink2)
  if (is.null(main)) main <- "Posterior medians, 50% and 95% intervals"
  par(xpd = NA)
  usr <- par("usr")
  left <- grconvertX(0.02, "ndc", "user")
  mtext(main, side = 3, line = 2, adj = 0, at = left, cex = 1.05, col = pal$ink)
  mtext(paste0("layer-specific parameters and hierarchical mean",
               if (!is.null(.scale_note(m))) paste0("  |  ", .scale_note(m)) else ""),
        side = 3, line = 0.8, adj = 0, at = left, cex = 0.7, col = pal$ink2)
  ents <- c(lapply(seq_len(K), function(k)
    list(type = "point", label = paste0("layer ", k, " (x", k - 1, " to x", k, ")"),
         pch = 21, bg = cols[k], pt.cex = 1.1)),
    list(list(type = "point", label = "hierarchical mean", pch = 23, bg = pal$ink, pt.cex = 1.2)))
  .legend_groups(list(list(title = "", entries = ents)), left,
                 usr[3] - 0.12 * diff(usr[3:4]), cex = 0.72)
  invisible(x)
}

plot.bergm_w <- function(x, what = c("phi", "mu", "sd"), type = c("trace", "density", "intervals"),
                         params = NULL, max_points = 2500, ...) {
  what <- match.arg(what)
  type <- match.arg(type)
  if (type == "intervals") return(.plot_intervals(x, ...))
  ml <- as.mcmc.list(x, what)
  labs <- colnames(ml[[1]])
  if (!is.null(params)) labs <- intersect(labs, params)
  if (!length(labs)) stop("no parameters to plot.", call. = FALSE)
  m <- x$model
  nc <- if (what == "phi" && is.null(params)) m$p else min(length(labs), 4L)
  nr <- ceiling(length(labs) / nc)
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  .ws_par(mfrow = c(nr, nc), mar = c(2.2, 3.2, 2, 0.6), oma = c(2.2, 0, 2.4, 0))
  cols <- .ws_pal$chains[(seq_along(ml) - 1L) %% length(.ws_pal$chains) + 1L]
  for (lab in labs) {
    vals <- lapply(ml, function(ch) as.numeric(ch[, lab]))
    if (type == "trace") {
      it <- as.numeric(stats::time(ml[[1]]))
      show <- if (length(it) > max_points) unique(round(seq(1, length(it), length.out = max_points))) else seq_along(it)
      plot(range(it), range(unlist(vals)), type = "n", xlab = "", ylab = "",
           axes = FALSE)
      abline(h = pretty(unlist(vals), 4), col = .ws_pal$grid, lwd = 0.6)
      for (c in seq_along(vals))
        lines(it[show], vals[[c]][show], col = adjustcolor(cols[c], 0.75), lwd = 0.6)
    } else {
      dens <- lapply(vals, stats::density)
      plot(range(sapply(dens, function(d) range(d$x))), c(0, max(sapply(dens, function(d) max(d$y)))),
           type = "n", xlab = "", ylab = "", axes = FALSE)
      abline(h = 0, col = .ws_pal$axis, lwd = 0.8)
      for (c in seq_along(dens)) lines(dens[[c]], col = cols[c], lwd = 1.4)
    }
    axis(1, col = .ws_pal$axis, lwd = 0, lwd.ticks = 0.6, cex.axis = 0.75)
    axis(2, col = .ws_pal$axis, lwd = 0, lwd.ticks = 0.6, cex.axis = 0.75)
    title(main = lab, cex.main = 0.85, adj = 0, line = 0.5)
  }
  mtext(if (type == "trace") "iteration" else "value", side = 1, outer = TRUE,
        line = 0.8, cex = 0.8, col = .ws_pal$ink2)
  mtext(paste0(if (type == "trace") "Traces" else "Posterior densities", " of ",
               switch(what, phi = "the layer-specific parameters", mu = "the hierarchical means",
                      sd = "the between-layer standard deviations")),
        side = 3, outer = TRUE, line = 0.9, adj = 0.01, cex = 0.85, col = .ws_pal$ink)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("topright", legend = paste("chain", seq_along(ml)), col = cols, lwd = 2,
         horiz = TRUE, bty = "n", cex = 0.75, text.col = .ws_pal$ink2, inset = c(0.01, 0.005))
  invisible(x)
}
