ws_prior <- function(model, mu0 = 0, Sigma0 = 3, nu0 = NULL, S0 = 1) {
  p <- if (inherits(model, "ws_model")) model$p else as.integer(model)
  terms <- if (inherits(model, "ws_model")) model$terms else paste0("s", seq_len(p))
  mu0 <- rep_len(as.numeric(mu0), p)
  as_mat <- function(M, name) {
    if (length(M) == 1) M <- diag(as.numeric(M), p)
    if (!is.matrix(M) || any(dim(M) != p))
      stop("'", name, "' must be a scalar or a ", p, " x ", p, " matrix.", call. = FALSE)
    M
  }
  Sigma0 <- as_mat(Sigma0, "Sigma0")
  S0 <- as_mat(S0, "S0")
  if (is.null(nu0)) nu0 <- p + 9
  if (nu0 <= p + 1) stop("'nu0' must exceed p + 1 = ", p + 1, ".", call. = FALSE)
  out <- list(mu0 = mu0, Sigma0 = Sigma0, nu0 = nu0, S0 = S0, p = p, terms = terms,
              Sigma0_inv = chol2inv(.chol_safe(Sigma0)),
              E_Sigma = S0 / (nu0 - p - 1))
  class(out) <- "ws_prior"
  out
}

print.ws_prior <- function(x, ...) {
  cat("Hierarchical prior (p = ", x$p, " statistics per layer)\n", sep = "")
  cat("  mu ~ N(mu0, Sigma0), Sigma ~ inverse-Wishart(nu0 = ", x$nu0, ", S0)\n", sep = "")
  sd_b <- sqrt(diag(x$E_Sigma))
  cat("  prior mean of mu0:", paste(.fmt(x$mu0, 2), collapse = ", "), "\n")
  cat("  prior sd of mu:", paste(.fmt(sqrt(diag(x$Sigma0)), 2), collapse = ", "), "\n")
  cat("  expected between-layer sd, sqrt(diag(E[Sigma])):",
      paste(.fmt(sd_b, 2), collapse = ", "), "\n")
  invisible(x)
}

.prior_cov_phi <- function(prior, K)
  kronecker(diag(K), prior$E_Sigma) + kronecker(matrix(1, K, K), prior$Sigma0)
