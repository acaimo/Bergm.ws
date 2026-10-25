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

plot.ws_predictive <- function(x, what = c("model", "gof"), style = c("intervals", "panels"),
                                main = NULL, ...) {
  what <- match.arg(what)
  style <- match.arg(style)
  if (what == "gof" && is.null(x$gof)) stop("no GOF statistics: rerun with gof = TRUE.", call. = FALSE)
  S <- if (what == "model") x$stats else x$gof
  obs <- if (what == "model") x$obs else x$gof_obs
  K <- x$model$K
  if (is.null(main))
    main <- paste0(if (x$type == "posterior") "Posterior" else "Prior", " predictive check: ",
                   if (what == "model") "model statistics" else "statistics not in the model")
  if (style == "panels") return(invisible(.pred_panels(S, obs, K, main)))
  invisible(.pred_intervals(S, obs, K, main))
}

.pred_intervals <- function(S, obs, K, main) {
  pal <- .ws_pal
  labs <- colnames(S)
  layer <- as.integer(sub("^L([0-9]+):.*$", "\\1", labs))
  stat <- sub("^L[0-9]+:", "", labs)
  m <- colMeans(S)
  s <- apply(S, 2, sd)
  s[s == 0] <- 1
  Z <- sweep(sweep(S, 2, m), 2, s, "/")
  z_obs <- (obs - m) / s
  q <- apply(Z, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
  pv <- .mid_p(S, obs)
  outside <- obs < apply(S, 2, quantile, 0.025) | obs > apply(S, 2, quantile, 0.975)

  y <- numeric(length(labs))
  pos <- 0
  for (k in seq_len(K)) {
    idx <- which(layer == k)
    pos <- pos + 1.2
    y[idx] <- -(pos + seq_along(idx) - 1)
    pos <- pos + length(idx)
  }
  xlim <- c(-4.2, 4.2)
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  .ws_par(mar = c(4, 9.5, 4.2, 5))
  plot.new()
  plot.window(xlim = xlim, ylim = c(min(y) - 0.8, 0.4))
  abline(v = c(-2, 2), col = pal$grid, lwd = 0.8)
  abline(v = 0, col = pal$axis, lwd = 1)
  clip <- function(v) pmin(pmax(v, xlim[1]), xlim[2])
  for (j in seq_along(labs)) {
    segments(clip(q[1, j]), y[j], clip(q[5, j]), y[j], col = pal$muted, lwd = 1.2)
    rect(clip(q[2, j]), y[j] - 0.18, clip(q[4, j]), y[j] + 0.18,
         col = "#b7d3f6", border = NA)
    segments(q[3, j], y[j] - 0.18, q[3, j], y[j] + 0.18, col = "#256abf", lwd = 1.5)
    zo <- clip(z_obs[j])
    points(zo, y[j], pch = if (outside[j]) 23 else 21, cex = 1.15, lwd = 1.5,
           bg = if (outside[j]) pal$critical else pal$ink, col = pal$surface)
    if (abs(z_obs[j]) > xlim[2])
      text(zo, y[j] + 0.38, sprintf("%.1f", z_obs[j]), cex = 0.6, col = pal$ink2)
  }
  axis(2, at = y, labels = stat, lwd = 0, cex.axis = 0.78, col.axis = pal$ink2)
  for (k in seq_len(K)) {
    idx <- which(layer == k)
    mtext(paste0("Layer ", k), side = 2, at = max(y[idx]) + 0.85, line = 8.8, adj = 0,
          cex = 0.8, col = pal$ink, font = 2)
  }
  axis(1, at = -4:4, lwd = 0, lwd.ticks = 0.6, cex.axis = 0.75)
  mtext("standardised value: (value - simulated mean) / simulated sd", side = 1,
        line = 2.3, cex = 0.75, col = pal$ink2)
  par(xpd = NA)
  usr <- par("usr")
  text(xlim[2] + 0.35, y, sprintf("%.2f", pv), adj = 0, cex = 0.65,
       col = ifelse(outside, pal$critical, pal$ink2))
  text(xlim[2] + 0.35, max(y) + 0.85, "mid-p", adj = 0, cex = 0.65, col = pal$muted)
  left <- grconvertX(0.02, "ndc", "user")
  mtext(main, side = 3, line = 2.7, adj = 0, at = left, cex = 0.95, col = pal$ink)
  .legend_groups(list(
    list(title = "simulated", entries = list(
      list(type = "box", label = "central 50% and median", col = "#b7d3f6", mid = "#256abf"),
      list(type = "line", label = "central 95%", col = pal$muted, lwd = 1.2, lty = "solid"))),
    list(title = "observed", entries = list(
      list(type = "point", label = "inside", pch = 21, bg = pal$ink, pt.cex = 1.1),
      list(type = "point", label = "outside the 95% interval", pch = 23, bg = pal$critical,
           pt.cex = 1.1)))),
    left, usr[4] + 0.045 * diff(usr[3:4]), cex = 0.7)
  invisible(data.frame(statistic = labs, z_observed = z_obs, mid_p = pv, outside = outside))
}

.pred_panels <- function(S, obs, K, main) {
  pal <- .ws_pal
  labs <- colnames(S)
  nc <- ncol(S) / K
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  .ws_par(mfrow = c(K, nc), mar = c(2, 1.2, 2, 0.6), oma = c(0.5, 0.5, 2.4, 0.5))
  for (j in seq_along(labs)) {
    v <- S[, j]
    br <- pretty(range(c(v, obs[j])), 20)
    h <- hist(v, breaks = br, plot = FALSE)
    plot(range(br), c(0, max(h$counts) * 1.12), type = "n", axes = FALSE, xlab = "", ylab = "")
    rect(h$breaks[-length(h$breaks)] + diff(h$breaks) * 0.06, 0,
         h$breaks[-1] - diff(h$breaks) * 0.06, h$counts, col = "#9ec5f4", border = NA)
    abline(h = 0, col = pal$axis, lwd = 0.8)
    abline(v = obs[j], col = pal$ink, lwd = 1.6)
    axis(1, lwd = 0, lwd.ticks = 0.6, cex.axis = 0.65)
    title(main = labs[j], cex.main = 0.75, adj = 0, line = 0.4, col.main = pal$ink2)
  }
  mtext(paste0(main, "  (bars: simulated, line: observed)"), side = 3, outer = TRUE,
        line = 0.8, adj = 0.01, cex = 0.85, col = pal$ink)
  invisible(NULL)
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
