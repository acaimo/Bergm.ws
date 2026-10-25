.chol_safe <- function(S, floor = 1e-10) {
  S <- (S + t(S)) / 2
  U <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(U)) {
    e <- eigen(S, symmetric = TRUE)
    v <- pmax(e$values, floor * max(1, max(abs(e$values))))
    S <- e$vectors %*% (v * t(e$vectors))
    U <- chol((S + t(S)) / 2)
  }
  U
}

.rmvn_u <- function(U) drop(rnorm(nrow(U)) %*% U)

.ldmvn_u <- function(x, mu, U) {
  z <- backsolve(U, x - mu, transpose = TRUE)
  -0.5 * sum(z^2) - sum(log(diag(U)))
}

.riwish <- function(nu, S) {
  W <- rWishart(1, nu, chol2inv(.chol_safe(S)))[, , 1]
  chol2inv(.chol_safe(W))
}

.vech <- function(M) M[lower.tri(M, diag = TRUE)]
.unvech <- function(v, p) {
  M <- matrix(0, p, p)
  M[lower.tri(M, diag = TRUE)] <- v
  M + t(M) - diag(diag(M), p)
}

.stats_vec <- function(res) drop(as.matrix(res$stats)[1, ])

.check_square_sym <- function(M, name) {
  if (!is.matrix(M) || nrow(M) != ncol(M))
    stop("'", name, "' must be a square matrix.", call. = FALSE)
  if (!isTRUE(all.equal(M, t(M), check.attributes = FALSE)))
    stop("'", name, "' must be symmetric (undirected networks only).", call. = FALSE)
  invisible(TRUE)
}

.fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
