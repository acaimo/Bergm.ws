ws_model <- function(data, formula, scale = NULL, labels = NULL) {
  if (!inherits(data, "ws_data")) stop("'data' must be created with ws_data().", call. = FALSE)
  if (missing(formula) || !inherits(formula, "formula"))
    stop("'formula' must be an ergm.multi formula on the layers x_1, ..., x_K",
         if (data$signed) " and the signs z", ".", call. = FALSE)
  K <- data$K
  rhs0 <- paste(deparse(formula[[length(formula)]], width.cutoff = 500L), collapse = " ")
  template <- grepl("\\bx_k\\b", rhs0, perl = TRUE)
  rhs <- if (template) {
    if (grepl("\\bx_[0-9]+\\b", rhs0, perl = TRUE))
      stop("a template formula uses x_k only, not numbered layers.", call. = FALSE)
    paste(vapply(seq_len(K), function(k) gsub("\\bx_k\\b", paste0("x_", k), rhs0, perl = TRUE), ""),
          collapse = " + ")
  } else rhs0
  if (!data$signed && grepl("\\bz\\b", rhs, perl = TRUE))
    stop("the formula uses the signs 'z', but the network has none.", call. = FALSE)
  f <- as.formula(paste("~", rhs), env = environment(formula))

  s_obs <- tryCatch(summary(f, basis = data$y), error = function(e)
    stop("the formula cannot be evaluated on the layers: ", conditionMessage(e),
         "\nLayers are named x_1, ..., x_", K, if (data$signed) " and the signs z",
         " (or x_k in a template formula).", call. = FALSE))
  nm <- names(s_obs)
  lay <- vapply(regmatches(nm, gregexpr("x_[0-9]+", nm)), function(u) {
    u <- unique(u)
    if (length(u) == 1) as.integer(sub("x_", "", u)) else NA_integer_
  }, integer(1))
  if (anyNA(lay) || any(!lay %in% seq_len(K)))
    stop("every statistic must use exactly one layer among x_1, ..., x_", K, "; not so for: ",
         paste(nm[is.na(lay) | !lay %in% seq_len(K)], collapse = ", "), call. = FALSE)
  if (is.unsorted(lay))
    stop("the statistics must be grouped by layer: all those of x_1, then those of x_2, ... ",
         "(a template formula with x_k does this automatically).", call. = FALSE)
  counts <- tabulate(lay, K)
  if (any(counts != counts[1]) || counts[1] == 0)
    stop("every layer needs the same number of statistics, in the same order; found ",
         paste(sprintf("x_%d: %d", seq_len(K), counts), collapse = ", "), ".", call. = FALSE)
  p <- counts[1]
  generic <- matrix(gsub("x_[0-9]+", "x_k", nm), K, p, byrow = TRUE)
  differ <- which(apply(generic, 2, function(v) length(unique(v)) > 1))
  if (length(differ))
    warning("statistic(s) ", paste(differ, collapse = ", "), " differ across layers (",
            paste(vapply(differ, function(j) paste(unique(generic[, j]), collapse = " vs "), ""),
                  collapse = "; "), "); the hierarchical prior pools them as if they were the ",
            "same statistic.", call. = FALSE)

  if (is.null(labels)) labels <- .short_labels(generic[1, ]) else {
    if (length(labels) != p || anyDuplicated(labels))
      stop("'labels' must be ", p, " distinct names, one per statistic of a layer.", call. = FALSE)
    labels <- as.character(labels)
  }
  sc <- .scale_vector(scale, nm, K, p)
  full_labels <- as.vector(t(outer(paste0("L", seq_len(K)), labels, paste, sep = ":")))
  names(s_obs) <- full_labels

  zero <- full_labels[s_obs == 0]
  if (length(zero))
    warning("statistic(s) equal to zero in the observed network: ", paste(zero, collapse = ", "),
            ". Their parameters are informed mainly by the prior.", call. = FALSE)

  out <- list(data = data, formula = f, rhs = rhs, template = if (template) rhs0 else NULL,
              terms = labels, labels = full_labels, K = K, p = p, signed = data$signed,
              scale = unname(sc), scale_type = if (is.character(scale)) scale else if (is.null(scale)) "none" else "user",
              s_obs = s_obs, s_obs_scaled = s_obs * sc, ergm_names = nm)
  class(out) <- "ws_model"

  mp <- .mple_design(out)
  if (!is.null(mp)) {
    out$change_cor <- lapply(seq_len(K), function(k) {
      Xk <- mp$predictor[, (k - 1L) * p + seq_len(p), drop = FALSE]
      keep <- rowSums(Xk != 0) > 0
      C <- .wcor(Xk[keep, , drop = FALSE], mp$weights[keep])
      dimnames(C) <- list(labels, labels)
      C
    })
    names(out$change_cor) <- paste0("layer ", seq_len(K))
    hi <- .high_cor(out$change_cor, 0.99)
    if (nrow(hi))
      warning("nearly collinear change statistics: ",
              paste(sprintf("%s: %s / %s (r = %.3f)", hi$layer, hi$a, hi$b, hi$r), collapse = "; "),
              ". Their parameters are separated mainly by the prior; consider dropping one.",
              call. = FALSE)
  }
  out
}

.scale_vector <- function(scale, nm, K, p) {
  if (is.null(scale)) return(rep(1, K * p))
  if (is.character(scale)) {
    if (!identical(scale, "gwesp")) stop("'scale' must be NULL, \"gwesp\" or numeric.", call. = FALSE)
    d <- suppressWarnings(as.numeric(sub("^.*~gwesp\\.fixed\\.([-0-9.eE+]+)$", "\\1", nm)))
    d[!grepl("~gwesp\\.fixed\\.", nm)] <- NA
    return(ifelse(is.na(d), 1, exp(-d)))
  }
  if (!length(scale) %in% c(p, K * p))
    stop("'scale' must have length ", p, " or ", K * p, ".", call. = FALSE)
  rep_len(as.numeric(scale), K * p)
}

.split_top <- function(s) {
  out <- character(0); depth <- 0; cur <- ""
  for (ch in strsplit(s, "")[[1]]) {
    if (ch == "," && depth == 0) { out <- c(out, cur); cur <- ""; next }
    if (ch == "(") depth <- depth + 1
    if (ch == ")") depth <- depth - 1
    cur <- paste0(cur, ch)
  }
  c(out, cur)
}

.logic_sign <- function(l) {
  l <- gsub("^\\((.*)\\)$", "\\1", gsub(" ", "", l))
  switch(l, "z&x_k" = "+", "!z&x_k" = "-", "x_k" = "", paste0("[", l, "]"))
}

.short_labels <- function(generic) {
  parts <- lapply(generic, function(g) {
    m <- regmatches(g, regexec("^L\\(pth=\\((.*)\\),bse=(.*),inord=(TRUE|FALSE)\\)~(.*)$", g))[[1]]
    if (length(m)) {
      paths <- vapply(.split_top(m[2]), .logic_sign, "")
      if (length(paths) == 1) paths <- rep(paths, 2)
      return(c(m[5], paste0(.logic_sign(m[3]), paste(paths, collapse = ""))))
    }
    m <- regmatches(g, regexec("^L\\((.*)\\)~(.*)$", g))[[1]]
    if (length(m)) return(c(m[3], .logic_sign(m[2])))
    c(g, "")
  })
  stat <- vapply(parts, `[`, "", 1)
  sgn <- vapply(parts, `[`, "", 2)
  fixed <- grepl("\\.fixed\\.[-0-9.eE+]+$", stat)
  base <- sub("\\.fixed\\.[-0-9.eE+]+$", "", stat)
  decay <- sub("^.*\\.fixed\\.", "", stat)
  short <- paste0(base, sgn)
  long <- ifelse(fixed, paste0(short, "(", decay, ")"), short)
  if (!anyDuplicated(short)) short else if (!anyDuplicated(long)) long else make.unique(long)
}

.wcor <- function(X, w) {
  w <- w / sum(w)
  Xc <- sweep(X, 2, colSums(X * w))
  C <- crossprod(Xc * w, Xc)
  s <- sqrt(pmax(diag(C), 0))
  R <- C / outer(s, s)
  R[!is.finite(R)] <- NA
  R
}

.high_cor <- function(cc, cut) {
  rows <- lapply(names(cc), function(nm) {
    C <- cc[[nm]]
    ij <- which(abs(C) > cut & upper.tri(C), arr.ind = TRUE)
    if (!nrow(ij)) return(NULL)
    data.frame(layer = nm, a = rownames(C)[ij[, 1]], b = colnames(C)[ij[, 2]],
               r = C[ij], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame(layer = character(0), a = character(0), b = character(0),
                               r = numeric(0)) else out
}

.scale_note <- function(m)
  switch(if (is.null(m$scale_type)) "none" else m$scale_type, gwesp = "GWESP statistics rescaled by exp(-decay)",
         user = "statistics multiplied by 'scale'", NULL)

print.ws_model <- function(x, ...) {
  cat("Multi-layer dissolution ERGM for", if (x$signed) "signed" else "unsigned",
      "weighted networks\n")
  cat("  layers (transitions):", x$K, "  statistics per layer:", x$p, "\n")
  cat("  formula", if (!is.null(x$template)) " (template, x_k = x_1, ..., x_K)", ":\n", sep = "")
  cat(strwrap(if (!is.null(x$template)) x$template else x$rhs, width = 88, indent = 4, exdent = 6),
      sep = "\n")
  if (!is.null(.scale_note(x))) cat(" ", .scale_note(x), "\n")
  cat("\nObserved statistics", if (any(x$scale != 1)) "(rescaled)", ":\n")
  print(round(matrix(x$s_obs_scaled, x$K, x$p, byrow = TRUE,
                     dimnames = list(paste0("layer ", seq_len(x$K)), x$terms)), 2))
  if (!is.null(x$change_cor)) {
    hi <- .high_cor(x$change_cor, 0.9)
    if (nrow(hi)) {
      cat("\nStrongly correlated change statistics (|r| > 0.9):\n")
      for (i in seq_len(nrow(hi)))
        cat(sprintf("  %s: %s / %s  r = %.3f\n", hi$layer[i], hi$a[i], hi$b[i], hi$r[i]))
    }
  }
  invisible(x)
}
