library(Bergm.ws)

sim <- ws_simulate(n = 24, seed = 5)
dat <- ws_data(sim$weights, sim$thresholds, sim$signs, sim$attributes)
stopifnot(inherits(dat, "ws_data"), dat$K == 3, dat$n == 24)
stopifnot(all(dat$levels[dat$levels >= 0] <= dat$K))
stopifnot(all(dat$signs[dat$levels >= 0] != 0), all(dat$signs[dat$levels < 0] == 0))

stopifnot(inherits(try(ws_data(sim$weights, rev(sim$thresholds), sim$signs), silent = TRUE),
                   "try-error"))

fs <- ~ L(~edges, ~(z & x_k)) +
  gwespL(0.7, fixed = TRUE, L.base = ~(z & x_k), Ls.path = list(~(z & x_k), ~(z & x_k))) +
  gwespL(0.7, fixed = TRUE, L.base = ~(!z & x_k), Ls.path = list(~(z & x_k), ~(z & x_k)))
mod <- suppressWarnings(ws_model(dat, fs, scale = "gwesp"))
stopifnot(mod$p == 3, length(mod$s_obs) == 9, identical(mod$terms, c("edges+", "gwesp+++", "gwesp-++")))
stopifnot(isTRUE(all.equal(unname(mod$scale), rep(c(1, exp(-0.7), exp(-0.7)), 3))))
stopifnot(inherits(try(ws_model(dat), silent = TRUE), "try-error"))
stopifnot(all(mapply(function(nm, k) grepl(paste0("x_", k), nm),
                     mod$ergm_names, rep(1:3, each = 3))))

pri <- ws_prior(mod)
stopifnot(pri$nu0 == mod$p + 9, isTRUE(all.equal(pri$E_Sigma, diag(3) / 8)))

U <- Bergm.ws:::.chol_safe(matrix(c(2, 1, 1, 2), 2))
stopifnot(isTRUE(all.equal(crossprod(U), matrix(c(2, 1, 1, 2), 2))))
v <- Bergm.ws:::.vech(matrix(c(2, 1, 1, 3), 2))
stopifnot(isTRUE(all.equal(Bergm.ws:::.unvech(v, 2), matrix(c(2, 1, 1, 3), 2))))

P <- matrix(0, 4, 4); N <- matrix(0, 4, 4)
P[1, 2] <- P[2, 3] <- P[1, 3] <- 1; P[2, 4] <- P[3, 4] <- 1; N[2, 3] <- 0
P <- P + t(P)
g <- Bergm.ws:::.gof_layer(P, N)
stopifnot(g[["triangles +++"]] == 2)
P2 <- P; P2[2, 3] <- P2[3, 2] <- 0; N2 <- N; N2[2, 3] <- N2[3, 2] <- 1
g2 <- Bergm.ws:::.gof_layer(P2, N2)
stopifnot(g2[["triangles +++"]] == 0, g2[["triangles -++"]] == 2)

set.seed(1)
fit <- bergm_ws(mod, n_iter = 120, burnin = 80, n_chains = 2, aux_iters = 500,
               n_info = 40, verbose = FALSE)
stopifnot(isTRUE(fit$settings$pool_adaptation))
stopifnot(isTRUE(all.equal(fit$chains[[1]]$gamma, fit$chains[[2]]$gamma)),
          isTRUE(all.equal(fit$chains[[1]]$shape, fit$chains[[2]]$shape)))
stopifnot(nrow(fit$chains[[1]]$accept) == 120, nrow(fit$chains[[1]]$phi) == 40)
s <- summary(fit)
stopifnot(nrow(s$phi) == 9, nrow(s$mu) == 3, all(is.finite(s$phi$mean)))
stopifnot(all(dim(coef(fit)) == c(3, 3)))

fitL <- bergm_ws(mod, update = "layer", n_iter = 30, burnin = 10, n_chains = 1,
                 aux_iters = 500, n_info = 40, collapse_mu = FALSE, delayed = FALSE,
                verbose = FALSE)
stopifnot(ncol(fitL$chains[[1]]$accept) == 3, all(fitL$chains[[1]]$n_aux == 1))

pri <- ws_prior(mod, Sigma0 = 2)
V <- kronecker(diag(3), pri$E_Sigma) + kronecker(matrix(1, 3, 3), pri$Sigma0)
set.seed(2)
draws <- t(replicate(4000, {
  mu <- drop(pri$mu0 + Bergm.ws:::.rmvn_u(chol(pri$Sigma0)))
  unlist(lapply(1:3, function(k) mu + Bergm.ws:::.rmvn_u(chol(pri$E_Sigma))))
}))
stopifnot(max(abs(cov(draws) - V)) < 0.25)
stopifnot(inherits(fit, "bergm_ws"))

u <- ws_simulate(n = 24, signed = FALSE, seed = 6)
stopifnot(is.null(u$signs))
datu <- ws_data(u$weights, u$thresholds, attributes = u$attributes)
stopifnot(!datu$signed, all(datu$signs[datu$levels >= 0] == 1))
fw <- ~ L(~edges, ~x_k) + L(~nodematch("module"), ~x_k) +
  gwespL(0.5, fixed = TRUE, L.base = ~x_k, Ls.path = list(~x_k, ~x_k)) +
  L(~gwdegree(0.5, fixed = TRUE), ~x_k)
modu <- suppressWarnings(ws_model(datu, fw))
stopifnot(modu$p == 4, identical(modu$terms, c("edges", "nodematch.module", "gwesp", "gwdeg")))
modl <- suppressWarnings(ws_model(datu, fw, labels = c("a", "b", "c", "d")))
stopifnot(identical(modl$terms, c("a", "b", "c", "d")))
stopifnot(inherits(try(ws_model(datu, fs), silent = TRUE), "try-error"))
stopifnot(inherits(try(bergm_ws(datu, fw), silent = TRUE), "try-error"))
fitu <- suppressWarnings(bergm_w(datu, fw, n_iter = 60, burnin = 20, n_chains = 2, aux_iters = 500,
                                 n_info = 40, verbose = FALSE))
stopifnot(identical(class(fitu), "bergm_w"), ncol(fitu$chains[[1]]$phi) == 12)
fitd <- suppressWarnings(bergm_w(dat, ~ L(~edges, ~x_k), n_iter = 30, burnin = 10, n_chains = 1, aux_iters = 500,
                                 n_info = 40, verbose = FALSE))
stopifnot(!fitd$model$signed)
pu <- ws_predictive(fitu, nsim = 5, burnin = 500, verbose = FALSE)
stopifnot(ncol(pu$gof) == 9)
f <- ~ L(~edges, ~(z & x_1)) + L(~triangle, ~(z & x_1)) +
       L(~edges, ~(z & x_2)) + L(~triangle, ~(z & x_2)) +
       L(~edges, ~(z & x_3)) + L(~triangle, ~(z & x_3))
mf <- ws_model(dat, formula = f)
stopifnot(mf$p == 2, is.null(mf$template), all(mf$scale == 1))
mt <- ws_model(dat, ~ L(~edges, ~(z & x_k)) + L(~triangle, ~(z & x_k)))
stopifnot(identical(unname(mt$s_obs), unname(mf$s_obs)), !is.null(mt$template))
mr <- suppressWarnings(ws_model(dat, as.formula(paste("~", mod$rhs)), scale = mod$scale[1:3]))
stopifnot(isTRUE(all.equal(unname(mr$s_obs_scaled), unname(mod$s_obs_scaled))))
stopifnot(inherits(try(ws_model(dat, ~ L(~edges, ~(z & x_k)) + L(~edges, ~(z & x_1))), silent = TRUE), "try-error"))
bad <- list(~ L(~edges, ~(z & x_2)) + L(~edges, ~(z & x_1)) + L(~edges, ~(z & x_3)),
            ~ L(~edges, ~(z & x_1)) + L(~edges, ~(z & x_2)),
            ~ L(~edges, ~(z & x_1 & x_2)) + L(~edges, ~(z & x_2)) + L(~edges, ~(z & x_3)))
for (b in bad) stopifnot(inherits(try(ws_model(dat, formula = b), silent = TRUE), "try-error"))
stopifnot(inherits(try(ws_model(datu, formula = f), silent = TRUE), "try-error"))
fu <- ~ L(~edges, ~x_1) + L(~edges, ~x_2) + L(~edges, ~x_3)
fitf <- suppressWarnings(bergm_w(datu, formula = fu, n_iter = 40, burnin = 10, n_chains = 1,
                                 aux_iters = 500, n_info = 40, verbose = FALSE))
stopifnot(ncol(fitf$chains[[1]]$phi) == 3)
cat("all tests passed\n")
