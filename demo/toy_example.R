library(Bergm.ws)

data(toy_ws)
dat <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)
dat

xy <- ws_layout(dat)
plot(dat, layout = xy, color_by = "module")
plot(dat, layout = "circle", color_by = "module")
plot_layers(dat, layout = xy, color_by = "module")
plot_adjacency(dat)

f <- ~ L(~edges, ~(z & x_k)) + L(~nodematch("pathway"), ~(z & x_k)) + L(~edges, ~(!z & x_k)) +
  gwespL(0.5, fixed = TRUE, L.base = ~(z & x_k), Ls.path = list(~(z & x_k), ~(z & x_k))) +
  gwespL(0.5, fixed = TRUE, L.base = ~(!z & x_k), Ls.path = list(~(z & x_k), ~(z & x_k)))
mod <- ws_model(dat, f, scale = "gwesp")
mod
pri <- ws_prior(mod)
plot(ws_predictive(mod, prior = pri, type = "prior", nsim = 100, seed = 1))

fit <- bergm_ws(mod, prior = pri, n_iter = 10000, burnin = 2500, thin = 2,
                n_chains = 4, cores = 2, seed = 2026)
summary(fit)
plot(fit, type = "intervals")
plot(fit)

post_pred <- ws_predictive(fit, nsim = 200, seed = 2)
post_pred
plot(post_pred, what = "gof")
ws_aux_check(fit, seed = 3)

datu <- ws_data(toy_ws$weights, toy_ws$thresholds, attributes = toy_ws$attributes)
plot(datu, color_by = "module")
fu <- ~ L(~edges, ~x_k) + L(~nodematch("pathway"), ~x_k) +
  gwespL(0.5, fixed = TRUE, L.base = ~x_k, Ls.path = list(~x_k, ~x_k))
fitu <- bergm_w(datu, fu, scale = "gwesp",
                n_iter = 6000, burnin = 1500, thin = 2, n_chains = 4, cores = 2, seed = 7)
summary(fitu)
plot(fitu, type = "intervals")
plot(ws_predictive(fitu, nsim = 200, seed = 8), what = "gof")
