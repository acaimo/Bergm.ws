library(Bergm.ws)
library(coda)
out_file <- "vignette_results.rds"
cores <- 2

data(toy_ws)

y <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)

s <- ~ L(~edges, ~(z & x_k)) + 
       L(~nodematch("pathway"), ~(z & x_k)) + 
       L(~edges, ~(!z & x_k)) +
       gwespL(0.5, fixed = TRUE, L.base = ~(z & x_k), Ls.path = list(~(z & x_k), ~(z & x_k))) +
       gwespL(0.5, fixed = TRUE, L.base = ~(!z & x_k), Ls.path = list(~(z & x_k), ~(z & x_k)))

model <- ws_model(y, s, scale = "gwesp")
prior <- ws_prior(model)

configs <- list(
  basic = list(delayed = FALSE, collapse_mu = FALSE, pool_adaptation = FALSE),
  ` + delayed acceptance` = list(delayed = TRUE, collapse_mu = FALSE, pool_adaptation = FALSE),
  ` + mu integrated out` = list(delayed = TRUE, collapse_mu = TRUE, pool_adaptation = FALSE),
  ` + pooled adaptation` = list(delayed = TRUE, collapse_mu = TRUE, pool_adaptation = TRUE))

rows <- list()
for (seed in c(11, 12)) for (nm in names(configs)) {
  post <- do.call(bergm_ws, c(list(model, 
                                   prior = prior, 
                                   n_iter = 5000, 
                                   burnin = 1000, 
                                   n_chains = 4,
                                   cores = cores, 
                                   seed = seed, 
                                   verbose = FALSE), configs[[nm]]))

  sp <- summary(post)
  ess <- c(sp$phi$ESS, sp$mu$ESS)
  
  rows[[length(rows) + 1]] <- data.frame(
    sampler = nm, seed = seed, seconds = post$elapsed, simulated = sp$sim_share,
    acceptance = mean(sp$acceptance), median_ESS = median(ess), min_ESS = min(ess),
    max_Rhat = max(c(sp$phi$Rhat, sp$mu$Rhat)), ESS_per_minute = 60 * median(ess) / post$elapsed)
}

comparison <- do.call(rbind, rows)

res <- list(comparison = comparison,
            info = list(R = R.version.string, 
                        date = Sys.Date(), 
                        cores = cores))

saveRDS(res, out_file, compress = "xz")
