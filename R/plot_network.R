ws_layout <- function(data, type = c("stress", "circle"), by = NULL,
                       iterations = 500, neg_length = 1.5) {
  if (!inherits(data, "ws_data")) stop("'data' must be created with ws_data().", call. = FALSE)
  type <- match.arg(type)
  X <- .stress_layout(data, iterations, neg_length)
  if (type == "circle") X <- .circle_layout(data, X, by)
  attr(X, "type") <- type
  X
}

.stress_layout <- function(data, iterations, neg_length) {
  L <- data$levels; S <- data$signs; n <- data$n
  D <- matrix(Inf, n, n)
  diag(D) <- 0
  e <- which(L >= 0, arr.ind = TRUE)
  D[e] <- ifelse(S[e] > 0, 1 / sqrt(1 + L[e]), neg_length)
  for (k in seq_len(n)) D <- pmin(D, outer(D[, k], D[k, ], "+"))
  if (any(!is.finite(D))) D[!is.finite(D)] <- 1.5 * max(D[is.finite(D)])
  X <- cmdscale(D, k = 2)
  if (ncol(X) < 2) X <- cbind(X, stats::rnorm(n, sd = 1e-3))
  W <- 1 / D^2
  diag(W) <- 0
  V <- -W
  diag(V) <- rowSums(W)
  J <- matrix(1 / n, n, n)
  Vp <- solve(V + J) - J
  for (t in seq_len(iterations)) {
    Dx <- as.matrix(stats::dist(X))
    Dx[Dx < 1e-9] <- 1e-9
    B <- -W * D / Dx
    diag(B) <- 0
    diag(B) <- -rowSums(B)
    Xn <- Vp %*% B %*% X
    done <- max(abs(Xn - X)) < 1e-7
    X <- Xn
    if (done) break
  }
  X <- sweep(X, 2, colMeans(X))
  X <- X / max(abs(X))
  dimnames(X) <- list(data$ids, c("x", "y"))
  X
}

.circle_layout <- function(data, X, by) {
  n <- data$n
  ang <- atan2(X[, 2], X[, 1])
  grp <- if (is.null(by)) rep("all", n) else {
    if (is.null(data$attributes) || !by %in% names(data$attributes))
      stop("attribute '", by, "' not found.", call. = FALSE)
    as.character(data$attributes[[by]])
  }
  g_ang <- tapply(seq_len(n), grp, function(i) atan2(mean(sin(ang[i])), mean(cos(ang[i]))))
  g_ord <- names(sort(g_ang, decreasing = TRUE))
  A <- pmax(data$levels + 1, 0) * (data$signs > 0)
  rank_in <- numeric(n)
  for (g in g_ord) {
    i <- which(grp == g)
    rel <- (ang[i] - g_ang[[g]] + pi) %% (2 * pi)
    if (length(i) > 2) {
      Ag <- A[i, i] + 1e-3
      Lg <- diag(rowSums(Ag)) - Ag
      f <- eigen(Lg, symmetric = TRUE)$vectors[, length(i) - 1]
      if (suppressWarnings(stats::cor(f, rel)) > 0) f <- -f
      rank_in[i] <- rank(f, ties.method = "first")
    } else rank_in[i] <- rank(-rel, ties.method = "first")
  }
  ord <- order(match(grp, g_ord), rank_in)
  g <- grp[ord]
  ng <- length(g_ord)
  gap <- if (ng > 1) 1.2 else 0
  pos <- cumsum(c(0, ifelse(g[-1] != g[-n], 1 + gap, 1)))
  total <- pos[n] + 1 + gap
  theta <- pi / 2 - 2 * pi * (pos + 0.5 * gap) / total
  Y <- matrix(NA_real_, n, 2, dimnames = list(data$ids, c("x", "y")))
  Y[ord, ] <- cbind(cos(theta), sin(theta))
  if (!is.null(by)) {
    attr(Y, "groups") <- grp
    attr(Y, "group_span") <- lapply(split(theta, g), range)
  }
  Y
}

.layout_of <- function(data, layout, by) {
  if (is.null(layout)) return(ws_layout(data))
  if (is.character(layout)) return(ws_layout(data, type = layout, by = by))
  xy <- as.matrix(layout)
  if (!is.numeric(xy) || nrow(xy) != data$n || ncol(xy) != 2)
    stop("'layout' must be \"stress\", \"circle\" or an n x 2 matrix.", call. = FALSE)
  if (is.null(attr(xy, "type"))) attr(xy, "type") <- "stress"
  xy
}

.bezier_q <- function(p0, p1, p2, npts = 18) {
  t <- seq(0, 1, length.out = npts)
  cbind((1 - t)^2 * p0[1] + 2 * (1 - t) * t * p1[1] + t^2 * p2[1],
        (1 - t)^2 * p0[2] + 2 * (1 - t) * t * p1[2] + t^2 * p2[2])
}

.edge_curve <- function(p0, p2, curvature, circle) {
  mid <- (p0 + p2) / 2
  if (circle) {
    chord <- sqrt(sum((p2 - p0)^2))
    return(.bezier_q(p0, mid * (0.75 - 0.35 * chord / 2), p2, npts = 28))
  }
  d <- p2 - p0
  len <- sqrt(sum(d^2))
  nrm <- if (len > 0) c(-d[2], d[1]) / len else c(0, 0)
  .bezier_q(p0, mid + curvature * len * nrm, p2)
}

.draw_edges <- function(levels, signs, xy, K, min_level = 0, curvature = 0.12,
                        neg_lty = "22", width = 1, uniform = FALSE, circle = FALSE) {
  e <- which(upper.tri(levels) & levels >= min_level, arr.ind = TRUE)
  if (!nrow(e)) return(invisible(NULL))
  lev <- levels[e]
  sg <- signs[e]
  for (l in sort(unique(lev))) for (s in c(-1, 1)) {
    sel <- which(lev == l & sg == s)
    if (!length(sel)) next
    pts <- do.call(rbind, lapply(sel, function(i)
      rbind(.edge_curve(xy[e[i, 1], ], xy[e[i, 2], ], curvature, circle), c(NA, NA))))
    frac <- if (uniform || K == 0) 0.6 else l / K
    col <- adjustcolor(if (s > 0) .ws_pal$pos else .ws_pal$neg,
                       alpha.f = 0.22 + 0.68 * frac)
    lines(pts, col = col, lwd = width * (0.55 + 2.4 * frac),
          lty = if (s < 0) neg_lty else 1, lend = 1)
  }
  invisible(NULL)
}

.node_style <- function(data, color_by) {
  n <- data$n
  if (is.null(color_by)) {
    return(list(bg = rep(.ws_pal$node, n), pch = rep(21, n), key = NULL, group_col = NULL))
  }
  if (is.null(data$attributes) || !color_by %in% names(data$attributes))
    stop("attribute '", color_by, "' not found.", call. = FALSE)
  v <- as.character(data$attributes[[color_by]])
  lev <- sort(unique(v))
  cols <- c(.ws_pal$node_cat, rep(.ws_pal$node_other, max(0, length(lev) - 3)))[seq_along(lev)]
  shapes <- c(21, 22, 24, 23)[(seq_along(lev) - 1) %% 4 + 1]
  key_lab <- if (length(lev) > 3) c(lev[1:3], "other") else lev
  key_col <- if (length(lev) > 3) c(cols[1:3], .ws_pal$node_other) else cols
  key_pch <- if (length(lev) > 3) c(shapes[1:3], 21) else shapes
  i <- match(v, lev)
  list(bg = cols[i], pch = ifelse(i > 3, 21, shapes[i]),
       key = list(label = key_lab, col = key_col, pch = key_pch, title = color_by),
       group_col = setNames(cols, lev))
}

.halo_text <- function(x, y, labels, cex = 0.62, col = .ws_pal$ink2, adj = c(0.5, 0.5),
                       srt = 0) {
  r <- strwidth("M", cex = cex) * 0.1
  for (a in seq(0, 2 * pi, length.out = 13)[-13])
    text(x + r * cos(a), y + r * sin(a), labels, adj = adj, srt = srt, cex = cex,
         col = .ws_pal$surface)
  text(x, y, labels, adj = adj, srt = srt, cex = cex, col = col)
}

.node_cex <- function(strength, size = 1)
  size * (0.85 + 1.5 * sqrt(strength / max(1, max(strength))))

.draw_nodes <- function(xy, strength, style, size = 1) {
  points(xy, pch = style$pch, bg = style$bg, col = .ws_pal$surface,
         cex = .node_cex(strength, size), lwd = 1.6)
}

.place_labels <- function(xy, idx, cex_nodes, cex = 0.62) {
  if (!length(idx)) return(matrix(numeric(0), 0, 2))
  u <- strheight("M", cex = 1)
  dirs <- seq(0, 2 * pi, length.out = 9)[-9]
  placed <- matrix(numeric(0), 0, 2)
  out <- matrix(NA_real_, length(idx), 2)
  for (j in seq_along(idx)) {
    i <- idx[j]
    r <- u * (0.35 * cex_nodes[i] + 0.7)
    cand <- cbind(xy[i, 1] + 1.25 * r * cos(dirs), xy[i, 2] + r * sin(dirs))
    others <- rbind(xy[-i, , drop = FALSE], placed)
    score <- apply(cand, 1, function(cc) min(sqrt((others[, 1] - cc[1])^2 + (others[, 2] - cc[2])^2)))
    score <- score + 0.05 * r * sin(dirs)
    out[j, ] <- cand[which.max(score), ]
    placed <- rbind(placed, out[j, ])
  }
  out
}

.which_labels <- function(labels, strength, circle = FALSE) {
  n <- length(strength)
  if (is.logical(labels)) return(rep_len(labels, n))
  labels <- match.arg(labels, c("auto", "none", "all"))
  switch(labels, none = rep(FALSE, n), all = rep(TRUE, n),
         auto = if (n <= 25 || (circle && n <= 80)) rep(TRUE, n) else
           rank(-strength, ties.method = "first") <= 6)
}

.circle_decor <- function(xy, lab, ids, style) {
  th <- atan2(xy[, 2], xy[, 1])
  span <- attr(xy, "group_span")
  r_lab <- 1.06
  if (!is.null(span) && !is.null(style$group_col)) {
    r_lab <- 1.1
    pad <- 0.42 * 2 * pi / length(th)
    for (g in names(span)) {
      tt <- seq(span[[g]][1] - pad, span[[g]][2] + pad, length.out = 60)
      col <- if (g %in% names(style$group_col)) style$group_col[[g]] else .ws_pal$node_other
      lines(1.045 * cos(tt), 1.045 * sin(tt), col = col, lwd = 3, lend = 1)
    }
  }
  for (i in which(lab)) {
    deg <- th[i] * 180 / pi
    left <- cos(th[i]) < -1e-8
    .halo_text(r_lab * cos(th[i]), r_lab * sin(th[i]), ids[i], cex = 0.55,
               adj = c(if (left) 1 else 0, 0.5), srt = if (left) deg + 180 else deg)
  }
}

.legend_groups <- function(groups, x0, y0, cex = 0.72, max_w = NULL) {
  pal <- .ws_pal
  if (is.null(max_w)) max_w <- par("usr")[2] - x0
  sw <- function(s) strwidth(s, cex = cex)
  em <- sw("M")
  key_w <- 2 * em
  row_h <- 2 * strheight("M", cex = cex)
  ent_w <- function(e) (if (e$type == "point") 0.7 * key_w else key_w) + 0.4 * em + sw(e$label) + em
  gw <- vapply(groups, function(g) sw(g$title) + 0.6 * em +
                 sum(vapply(g$entries, ent_w, numeric(1))), numeric(1))
  x <- x0; y <- y0
  for (j in seq_along(groups)) {
    if (x > x0 && x + gw[j] > x0 + max_w) { x <- x0; y <- y - row_h }
    g <- groups[[j]]
    text(x, y, g$title, adj = c(0, 0.5), cex = cex, col = pal$muted)
    x <- x + sw(g$title) + 0.6 * em
    for (e in g$entries) {
      if (e$type == "line") {
        segments(x, y, x + key_w, y, col = e$col, lwd = e$lwd, lty = e$lty, lend = 1)
        x <- x + key_w + 0.4 * em
      } else if (e$type == "box") {
        hh <- 0.45 * strheight("M", cex = cex)
        rect(x, y - hh, x + key_w, y + hh, col = e$col, border = NA)
        if (!is.null(e$mid))
          segments(x + key_w / 2, y - hh, x + key_w / 2, y + hh, col = e$mid, lwd = 1.5)
        x <- x + key_w + 0.4 * em
      } else {
        points(x + 0.35 * key_w, y, pch = e$pch, bg = e$bg,
               col = if (is.null(e$col)) pal$surface else e$col, cex = e$pt.cex, lwd = 1.2)
        x <- x + 0.7 * key_w + 0.4 * em
      }
      text(x, y, e$label, adj = c(0, 0.5), cex = cex, col = pal$ink2)
      x <- x + sw(e$label) + em
    }
    x <- x + 1.4 * em
  }
  invisible(y)
}

.network_legend <- function(K, neg_lty, style, signed = TRUE) {
  pal <- .ws_pal
  fr <- (0:K) / max(1, K)
  groups <- list(
    list(title = "sign", entries = list(
      list(type = "line", label = "positive", col = pal$pos, lwd = 2.2, lty = "solid"),
      list(type = "line", label = "negative", col = pal$neg, lwd = 2.2, lty = neg_lty))),
    list(title = "weight level", entries = lapply(0:K, function(l)
      list(type = "line", label = paste0("x", l), lwd = 0.55 + 2.4 * fr[l + 1], lty = "solid",
           col = adjustcolor(pal$ink2, alpha.f = 0.3 + 0.6 * fr[l + 1])))))
  if (!signed) groups[[1]] <- NULL
  if (!is.null(style$key))
    groups[[length(groups) + 1]] <- list(title = style$key$title, entries = lapply(seq_along(style$key$label), function(i)
      list(type = "point", label = style$key$label[i], pch = style$key$pch[i],
           bg = style$key$col[i], pt.cex = 1.3)))
  par(xpd = NA)
  usr <- par("usr")
  .legend_groups(groups, usr[1], usr[3] - 0.035 * diff(usr[3:4]))
}

plot.ws_data <- function(x, layout = NULL, color_by = NULL,
                          labels = c("auto", "none", "all"), main = NULL,
                          curvature = 0.12, neg_lty = "22", legend = TRUE, ...) {
  pal <- .ws_pal
  K <- x$K
  xy <- .layout_of(x, layout, color_by)
  circle <- identical(attr(xy, "type"), "circle")
  strength <- rowSums(pmax(x$levels + 1, 0))
  style <- .node_style(x, color_by)
  lab <- .which_labels(labels, strength, circle)
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  .ws_par(mar = c(if (legend) 4.6 else 0.6, 0.6, 3.6, 0.6))
  plot.new()
  pad <- if (circle) 0.32 else 0.1
  plot.window(xlim = range(xy[, 1]) + c(-pad, pad), ylim = range(xy[, 2]) + c(-pad, pad), asp = 1)
  .draw_edges(x$levels, x$signs, xy, K, 0, curvature, neg_lty, circle = circle)
  cex_nodes <- .node_cex(strength)
  .draw_nodes(xy, strength, style)
  if (circle) {
    .circle_decor(xy, lab, x$ids, style)
  } else if (any(lab)) {
    idx <- which(lab)
    idx <- idx[order(-strength[idx])]
    lp <- .place_labels(xy, idx, cex_nodes)
    .halo_text(lp[, 1], lp[, 2], x$ids[idx])
  }

  E <- sum(x$levels >= 0) / 2
  pos <- sum(x$levels >= 0 & x$signs > 0) / 2
  if (is.null(main)) main <- if (x$signed) "Signed weighted network" else "Weighted network"
  title(main = main, adj = 0, line = 2.2, cex.main = 1.15)
  mtext(sprintf("%d nodes  |  %d edges in x0%s  |  node size: weighted degree", x$n, E,
                if (x$signed) sprintf(", %.0f%% positive", 100 * pos / max(1, E)) else ""),
        side = 3, line = 1.0, adj = 0, cex = 0.72, col = pal$ink2)
  if (legend) .network_legend(K, neg_lty, style, x$signed)
  invisible(xy)
}

plot_layers <- function(data, layout = NULL, color_by = NULL, ncol = NULL,
                        curvature = 0.12, neg_lty = "22", main = NULL) {
  if (!inherits(data, "ws_data")) stop("'data' must be created with ws_data().", call. = FALSE)
  pal <- .ws_pal
  K <- data$K
  xy <- .layout_of(data, layout, color_by)
  circle <- identical(attr(xy, "type"), "circle")
  style <- .node_style(data, color_by)
  if (is.null(ncol)) ncol <- K + 1
  nrow <- ceiling((K + 1) / ncol)
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  .ws_par(mfrow = c(nrow, ncol), mar = c(0.4, 0.4, 3.2, 0.4), oma = c(0, 0, 2.2, 0))
  pad <- if (circle) 0.12 else 0.1
  for (k in 0:K) {
    plot.new()
    plot.window(xlim = range(xy[, 1]) + c(-pad, pad), ylim = range(xy[, 2]) + c(-pad, pad), asp = 1)
    .draw_edges(data$levels, data$signs, xy, K, min_level = k, curvature, neg_lty,
                uniform = TRUE, circle = circle)
    active <- rowSums(data$levels >= k) > 0
    st <- style
    st$bg[!active] <- "#dcdad3"
    deg <- rowSums(data$levels >= k)
    .draw_nodes(xy, deg, st, size = 0.8)
    E <- sum(data$levels >= k) / 2
    pos <- sum(data$levels >= k & data$signs > 0) / 2
    title(main = bquote(x[.(k)]), adj = 0, line = 1.6, cex.main = 1.1)
    mtext(sprintf("weight >= %s  |  %d edges%s", format(data$thresholds[k + 1], digits = 3), E,
                  if (data$signed) sprintf(", %.0f%% positive", 100 * pos / max(1, E)) else ""),
          side = 3, line = 0.5, adj = 0, cex = 0.62, col = pal$ink2)
  }
  if (is.null(main)) main <- if (data$signed)
    "Nested layers: positive (blue) and negative (red, dashed) edges" else "Nested layers"
  mtext(main, side = 3, outer = TRUE, line = 0.6, adj = 0.01, cex = 0.85, col = pal$ink)
  invisible(xy)
}

plot_adjacency <- function(data, order = c("layout", "attribute", "none"), by = NULL,
                           labels = NULL, main = NULL) {
  if (!inherits(data, "ws_data")) stop("'data' must be created with ws_data().", call. = FALSE)
  order <- match.arg(order)
  pal <- .ws_pal
  K <- data$K; n <- data$n
  V <- ifelse(data$levels >= 0, data$signs * (data$levels + 1), 0)
  strength <- rowSums(pmax(data$levels + 1, 0))
  ord <- switch(order,
    none = seq_len(n),
    layout = {
      xy <- ws_layout(data)
      hclust(stats::dist(xy), method = "average")$order
    },
    attribute = {
      if (is.null(by) || is.null(data$attributes) || !by %in% names(data$attributes))
        stop("order = \"attribute\" needs 'by', the name of a node attribute.", call. = FALSE)
      order(as.character(data$attributes[[by]]), -strength)
    })
  M <- V[ord, ord]
  L <- K + 1
  cols <- c(rev(.ws_level_cols(L, "neg")), pal$neutral, .ws_level_cols(L, "pos"))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  if (is.null(labels)) labels <- n <= 50
  .ws_par(mar = c(if (labels) 3 else 1, if (labels) 3 else 1, 3.4, 7.5))
  image(seq_len(n), seq_len(n), M[, n:1], col = cols, breaks = seq(-L - 0.5, L + 0.5, by = 1),
        axes = FALSE, xlab = "", ylab = "", asp = 1, useRaster = TRUE)
  if (labels) {
    par(xpd = NA)
    off <- 0.25 + 0.3 * n / 40
    text(seq_len(n), 0.5 - off, data$ids[ord], srt = 90, adj = c(1, 0.5), cex = 0.45,
         col = pal$muted)
    text(0.5 - off, seq_len(n), rev(data$ids[ord]), adj = c(1, 0.5), cex = 0.45,
         col = pal$muted)
    par(xpd = FALSE)
  }
  if (is.null(main)) main <- if (data$signed) "Signed adjacency by weight level" else "Adjacency by weight level"
  title(main = main, adj = 0, line = 1.8, cex.main = 1.05)
  mtext("darker = edge still present in higher layers (larger weight)", side = 3, line = 0.6,
        adj = 0, cex = 0.7, col = pal$ink2)
  par(xpd = NA)
  usr <- par("usr")
  x0 <- usr[2] + 0.04 * diff(usr[1:2])
  h <- diff(usr[3:4]) / (2 * L + 6)
  yy <- usr[4] - h * 0.5
  text(x0, yy, if (data$signed) "positive" else "edges", adj = c(0, 0.5), cex = 0.7, col = pal$ink2)
  for (l in L:1) {
    yy <- yy - h
    rect(x0, yy - h * 0.4, x0 + h * 0.9, yy + h * 0.4, col = .ws_level_cols(L, "pos")[l], border = NA)
    text(x0 + h * 1.2, yy, .level_label(l), adj = c(0, 0.5), cex = 0.62, col = pal$ink2)
  }
  if (!data$signed) return(invisible(ord))
  yy <- yy - 1.4 * h
  text(x0, yy, "negative", adj = c(0, 0.5), cex = 0.7, col = pal$ink2)
  for (l in L:1) {
    yy <- yy - h
    rect(x0, yy - h * 0.4, x0 + h * 0.9, yy + h * 0.4, col = .ws_level_cols(L, "neg")[l], border = NA)
    text(x0 + h * 1.2, yy, .level_label(l), adj = c(0, 0.5), cex = 0.62, col = pal$ink2)
  }
  invisible(ord)
}

.level_label <- function(l) if (l == 1) "x0 only" else paste0("x0 to x", l - 1)
