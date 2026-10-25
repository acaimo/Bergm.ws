# Bergm.ws 0.1.0

* First version.
* `bergm_w()` fits multi-layer dissolution ERGMs to weighted networks;
  `bergm_ws()` to signed weighted networks, with the signs modelled jointly
  with the layers.
* Data: nested layers from weights and thresholds, with optional signs
  (`ws_data()`); a simulator of small weighted networks, signed or not
  (`ws_simulate()`), and a simulated example (`toy_ws`).
* Model: an `ergm.multi` formula on the layers (`ws_model()`), with checks that every statistic is computed on its own layer, that no
  statistic is zero in the data and that no two statistics have nearly
  collinear change statistics.
* Normal-inverse-Wishart hierarchical prior (`ws_prior()`).
* Adaptive exchange algorithm: joint or layer-wise updates, penalised MPLE
  start with moment-matching Newton steps, delayed acceptance, hierarchical
  mean integrated out of the exchange step, proposal covariance learned during
  the burn-in, pooled across chains and frozen, chains in parallel.
* Summaries with Rhat and effective sample sizes; trace, density and interval
  plots.
* Prior and posterior predictive checks with statistics not in the model, and
  a check of the length of the auxiliary chains.
* Network plots: stress and circular layouts, layers side by side, adjacency
  matrix by weight level.
* One Quarto vignette (`Bergm.ws.qmd`): model and notation, signed networks
  with `bergm_ws()`, networks without signs with `bergm_w()`, and the adaptive
  exchange algorithm; the code uses the names of the notation.
* Models are written only as `ergm.multi` formulas (`ws_model()`, `bergm_w()`,
  `bergm_ws()`), layer by layer or as a template with `x_k`; they are checked
  layer by layer, statistics get short labels, and `scale = "gwesp"` rescales
  GWESP statistics.
