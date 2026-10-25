ws_simulate <- function(n = 40, K = 3, n_modules = 3, density = 0.30,
                        retention = c(0.55, 0.50, 0.50), signed = TRUE,
                        faction_share = 0.18, sign_noise = 0.08, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n >= 10, K >= 1, n_modules >= 1, density > 0, density < 1)
  retention <- rep_len(retention, K)

  module <- sample(rep_len(seq_len(n_modules), n))
  ang <- 2 * pi * (seq_len(n_modules) - 1) / n_modules
  centre <- cbind(cos(ang), sin(ang)) * 2.2
  pos <- centre[module, , drop = FALSE] + matrix(rnorm(2 * n, sd = 0.75), n, 2)

  faction <- ifelse(runif(n) < faction_share, -1, 1)
  between <- outer(faction, faction) < 0

  ut <- upper.tri(diag(n))
  noise <- matrix(0, n, n)
  noise[ut] <- rnorm(sum(ut), sd = 0.25)
  noise <- noise + t(noise)
  score <- -as.matrix(stats::dist(pos)) + noise - if (signed) 0.8 * between else 0
  diag(score) <- -Inf

  w <- matrix(0, n, n)
  w[ut] <- rank(score[ut]) / (sum(ut) + 1)
  w <- w + t(w)
  keep <- density * cumprod(c(1, retention))
  thresholds <- as.numeric(stats::quantile(w[ut], 1 - keep, names = FALSE))

  sgn <- NULL
  if (signed) {
    sgn <- outer(faction, faction)
    p_flip <- sign_noise * 2 * (1 - pmin(1, pmax(0, (w - thresholds[1]) / (1 - thresholds[1]))))
    flip <- matrix(0, n, n)
    flip[ut] <- rbinom(sum(ut), 1, p_flip[ut])
    flip <- flip + t(flip)
    sgn[flip == 1] <- -sgn[flip == 1]
    diag(sgn) <- 0
  }
  ids <- sprintf("g%02d", seq_len(n))
  dimnames(w) <- list(ids, ids)
  if (signed) dimnames(sgn) <- list(ids, ids)
  p_path <- seq(0.3, 0.7, length.out = n_modules)[module]
  attributes <- data.frame(module = paste0("M", module),
                           pathway = ifelse(runif(n) < p_path, "yes", "no"),
                           row.names = ids, stringsAsFactors = FALSE)
  out <- list(weights = round(w, 4), signs = sgn, thresholds = round(thresholds, 4),
              attributes = attributes)
  if (!signed) out$signs <- NULL
  out
}

.levels_from_weights <- function(weights, thresholds) {
  L <- matrix(-1L, nrow(weights), ncol(weights))
  for (k in seq_along(thresholds)) L[weights >= thresholds[k]] <- k - 1L
  diag(L) <- -1L
  L[is.na(weights)] <- -1L
  L
}

.make_net <- function(adj, attributes, ids) {
  g <- network::network(adj, directed = FALSE, matrix.type = "adjacency")
  network::set.vertex.attribute(g, "vertex.names", ids)
  if (!is.null(attributes))
    for (a in names(attributes))
      network::set.vertex.attribute(g, a, attributes[[a]])
  g
}

ws_data <- function(weights, thresholds, signs = NULL, attributes = NULL) {
  .check_square_sym(weights, "weights")
  signed <- !is.null(signs)
  if (signed) {
    .check_square_sym(signs, "signs")
    if (!all(dim(weights) == dim(signs)))
      stop("'weights' and 'signs' must have the same dimensions.", call. = FALSE)
  }
  thresholds <- as.numeric(thresholds)
  if (length(thresholds) < 2 || is.unsorted(thresholds, strictly = TRUE))
    stop("'thresholds' must be strictly increasing, with at least two values ",
         "(baseline layer x0 and one higher layer).", call. = FALSE)
  n <- nrow(weights)
  K <- length(thresholds) - 1L
  ids <- rownames(weights)
  if (is.null(ids)) ids <- sprintf("v%02d", seq_len(n))
  if (!is.null(attributes)) {
    attributes <- as.data.frame(attributes, stringsAsFactors = FALSE)
    if (nrow(attributes) != n) stop("'attributes' must have one row per node.", call. = FALSE)
  }

  levels <- .levels_from_weights(weights, thresholds)
  s <- if (signed) sign(signs) else (levels >= 0) * 1
  s[levels < 0] <- 0
  if (signed && any(levels >= 0 & s == 0))
    stop("every edge of the baseline layer needs a non-zero sign.", call. = FALSE)
  dimnames(levels) <- dimnames(s) <- list(ids, ids)

  x <- lapply(0:K, function(k) .make_net((levels >= k) * 1, attributes, ids))
  names(x) <- paste0("x_", 0:K)
  z <- if (signed) .make_net(((levels >= 0) & s > 0) * 1, attributes, ids) else NULL
  y <- Layer(if (signed) c(list(z = z), x[-1]) else x[-1])

  env <- new.env(parent = asNamespace("Bergm.ws"))
  for (k in 0:(K - 1)) assign(paste0("x_", k), x[[k + 1]], envir = env)
  supports <- c(if (signed) "x_0", paste0("x_", 0:(K - 1)))
  constraint <- as.formula(sprintf("~Dyads(vary = ~dyadcov(Layer(list(%s))))",
                                   paste(supports, collapse = ", ")), env = env)

  out <- list(y = y, x = x, z = z, K = K, n = n, ids = ids, signed = signed,
              levels = levels, signs = s, weights = weights, thresholds = thresholds,
              attributes = attributes, constraint = constraint)
  class(out) <- "ws_data"
  out
}

.layer_adj <- function(levels, signs, k) {
  A <- (levels >= k) * 1
  list(P = A * (signs > 0), N = A * (signs < 0))
}

.gof_layer <- function(P, N, signed = TRUE) {
  P2 <- P %*% P
  if (!signed)
    return(c(triangles = sum(diag(P2 %*% P)) / 6,
             `sd(degree)` = stats::sd(rowSums(P)),
             `max(degree)` = max(rowSums(P))))
  N2 <- N %*% N
  c(`triangles +++` = sum(diag(P2 %*% P)) / 6,
    `triangles -++` = sum(N * P2) / 2,
    `triangles +--` = sum(P * N2) / 2,
    `triangles ---` = sum(diag(N2 %*% N)) / 6,
    `sd(degree+)`   = stats::sd(rowSums(P)))
}

.gof_stats <- function(signed)
  if (signed) c("triangles +++", "triangles -++", "triangles +--", "triangles ---", "sd(degree+)") else
    c("triangles", "sd(degree)", "max(degree)")

.gof_names <- function(K, signed = TRUE)
  as.vector(t(outer(paste0("L", seq_len(K)), .gof_stats(signed), paste, sep = ":")))

.gof_observed <- function(data) {
  g <- unlist(lapply(seq_len(data$K), function(k) {
    a <- .layer_adj(data$levels, data$signs, k)
    .gof_layer(a$P, a$N, data$signed)
  }))
  setNames(g, .gof_names(data$K, data$signed))
}

.gof_simulated <- function(nw, K, signed = TRUE) {
  layers <- uncombine_network(nw)
  adj <- function(g) as.matrix(g, matrix.type = "adjacency")
  zmat <- if (signed) adj(layers[["z"]]) else NULL
  g <- unlist(lapply(seq_len(K), function(k) {
    A <- adj(layers[[paste0("x_", k)]])
    if (signed) .gof_layer(A * zmat, A * (1 - zmat)) else .gof_layer(A, 0 * A, FALSE)
  }))
  setNames(g, .gof_names(K, signed))
}

.layer_table <- function(data) {
  rows <- lapply(0:data$K, function(k) {
    a <- .layer_adj(data$levels, data$signs, k)
    e_pos <- sum(a$P) / 2
    e_neg <- sum(a$N) / 2
    g <- .gof_layer(a$P, a$N, data$signed)
    if (!data$signed)
      return(data.frame(layer = paste0("x_", k), threshold = data$thresholds[k + 1],
                        edges = e_pos, triangles = g[["triangles"]],
                        `mean degree` = round(2 * e_pos / data$n, 2), check.names = FALSE))
    tri <- sum(g[1:4])
    data.frame(layer = paste0("x_", k), threshold = data$thresholds[k + 1],
               edges = e_pos + e_neg, positive = e_pos, negative = e_neg,
               `% positive` = round(100 * e_pos / max(1, e_pos + e_neg), 1),
               triangles = tri, `% +++` = round(100 * g[1] / max(1, tri), 1),
               check.names = FALSE)
  })
  do.call(rbind, rows)
}

print.ws_data <- function(x, ...) {
  cat(if (x$signed) "Signed" else "Unsigned", " weighted network: ", x$n, " nodes, K = ", x$K,
      " layer transitions (x_0 -> x_", x$K, ")\n", sep = "")
  if (!is.null(x$attributes))
    cat("Node attributes:", paste(names(x$attributes), collapse = ", "), "\n")
  cat("\n")
  print(.layer_table(x), row.names = FALSE)
  invisible(x)
}
