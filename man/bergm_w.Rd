\name{bergm_w}
\alias{bergm_w}
\alias{bergm_ws}
\title{Bayesian multi-layer dissolution ERGMs for weighted and signed weighted networks}
\description{
\code{bergm_w()} fits the model to a weighted network without signs, with
statistics computed on each layer \eqn{x_k}; \code{bergm_ws()} fits it to a
signed weighted network, with the signs modelled jointly with the layers and
statistics computed on the positive and negative edges of each layer.

Both sample the posterior of the layer-specific parameters \eqn{\varphi_k} and of
the hierarchical parameters \eqn{(\mu, \Sigma)} with an adaptive approximate
exchange algorithm. Each exchange step simulates an auxiliary network from the
joint model on the observed supports, starting at the observed network, so
that the intractable normalising constant cancels; \eqn{\mu} and \eqn{\Sigma}
are then drawn from their Normal and inverse-Wishart full conditionals.

Compared with the basic algorithm, proposals are screened by delayed
acceptance, the hierarchical mean is integrated out of the exchange step, and
the proposal covariance is learned from the chains during the burn-in, pooled
across chains and then frozen. None of these changes the posterior.
}
\usage{
bergm_w(data, formula = NULL, scale = NULL, labels = NULL, prior = NULL,
        update = c("joint", "layer"), n_iter = 10000,
        burnin = floor(n_iter / 4), thin = 1, n_chains = 4, cores = 1,
        aux_iters = NULL, init = "mple", prop_cov = "information",
        adapt_cov = TRUE, pool_adaptation = TRUE, adapt_every = NULL,
        adapt_start = NULL, adapt_scaling = 0.1, target_accept = 0.234,
        delayed = TRUE, collapse_mu = TRUE, n_info = 200, n_newton = 2,
        dispersion = 4, seed = NULL, verbose = TRUE)

bergm_ws(data, formula = NULL, scale = NULL, labels = NULL, prior = NULL,
         update = c("joint", "layer"), n_iter = 10000,
         burnin = floor(n_iter / 4), thin = 1, n_chains = 4, cores = 1,
         aux_iters = NULL, init = "mple", prop_cov = "information",
         adapt_cov = TRUE, pool_adaptation = TRUE, adapt_every = NULL,
         adapt_start = NULL, adapt_scaling = 0.1, target_accept = 0.234,
         delayed = TRUE, collapse_mu = TRUE, n_info = 200, n_newton = 2,
         dispersion = 4, seed = NULL, verbose = TRUE)
}
\arguments{
  \item{data}{a \code{\link{ws_data}} object, or a \code{\link{ws_model}}
    (then \code{formula}, \code{scale} and \code{labels} are not needed).
    \code{bergm_w()} uses only the weights of signed data;
    \code{bergm_ws()} needs signed data.}
  \item{formula, scale, labels}{the \pkg{ergm.multi} formula on the layers
    \code{x_1}, \dots, \code{x_K} (or the template with \code{x_k}) and, for
    \code{bergm_ws()}, the signs \code{z}; multipliers and names of the
    statistics. See \code{\link{ws_model}}.}
  \item{prior}{a \code{\link{ws_prior}}; default \code{ws_prior(model)}.}
  \item{update}{\code{"joint"}: one exchange step for all \eqn{Kp}
    parameters per iteration, with one auxiliary network; robust to posterior
    correlations both within and between layers. \code{"layer"}: one exchange
    step per layer with the other layers held fixed (\eqn{K} auxiliary
    networks per iteration); adequate when parameters of different layers are
    nearly uncorrelated.}
  \item{n_iter, burnin, thin}{iterations per chain, burn-in and thinning.}
  \item{n_chains}{number of chains. The first starts at the initial value,
    the others at draws from \eqn{N(\varphi^{(0)}, \code{dispersion} \times
    B^{(0)})}, so that convergence can be checked with \eqn{\hat R}.}
  \item{cores}{chains run in parallel with \code{parallel::mclapply} when
    \code{cores > 1} (not on Windows).}
  \item{aux_iters}{MCMC steps for each auxiliary network. The default is ten
    proposals per dyad-layer pair that can change (at least 2000); check it
    with \code{\link{ws_aux_check}}.}
  \item{init}{\code{"mple"} (default): penalised maximum pseudo-likelihood,
    i.e. the logistic pseudo-likelihood times the Gaussian prior, which is
    always finite; \code{"zero"}; or a numeric vector.}
  \item{prop_cov}{initial proposal covariance \eqn{B^{(0)}}:
    \code{"information"} (default) uses
    \eqn{\{\mathrm{Cov}(s) + P\}^{-1}}, with \eqn{\mathrm{Cov}(s)} the
    covariance of statistics simulated at the starting value (the Fisher
    information) and \eqn{P} the prior precision; \code{"diagonal"} uses
    \eqn{0.01 I}; or a matrix.}
  \item{adapt_cov}{if \code{TRUE}, the proposal covariance is re-estimated
    from the second half of the chain's history at every adaptation step
    (adaptive Metropolis), mixed 95:5 with \eqn{B^{(0)}} so that it cannot
    collapse.}
  \item{pool_adaptation}{with several chains and \code{adapt_cov = TRUE}, the
    burn-in has three phases. (1) During the first three quarters every chain
    adapts its own proposal covariance and scale. (2) The within-chain
    covariances are averaged, the pooled matrix is given to every chain, and
    during the last quarter of the burn-in the chains only tune the proposal
    scale. (3) The scale is set to the geometric mean of the chains' scales and
    the proposal is frozen. The retained draws of all chains then come from the
    same Markov kernel, with a covariance estimated from four times as many
    draws, and a chain that adapted to a poor local shape early on is
    corrected.}
  \item{adapt_every, adapt_start}{adaptation interval and first adaptation
    iteration; by default \code{max(20, burnin \%/\% 20)} and twice that.
    Adaptation stops at the end of the burn-in, so the retained draws come from
    a Markov chain with a fixed kernel.}
  \item{adapt_scaling}{the proposal scale is multiplied by
    \eqn{1 \pm} \code{adapt_scaling} according to the acceptance rate in the
    last window.}
  \item{target_accept}{target acceptance rate.}
  \item{delayed}{delayed acceptance (Christen and Fox, 2005). Each proposal is
    first screened with a Gaussian approximation of the log-likelihood,
    \eqn{g^\top(\varphi - c) - (\varphi - c)^\top H (\varphi - c)/2}, where
    \eqn{H} and \eqn{g} are the Fisher information and the score at the
    starting value \eqn{c}, together with the exact prior; the auxiliary
    network is simulated only for proposals that pass, and a second stage
    with the exchange ratio divided by the approximate likelihood ratio keeps
    the posterior exact. Most proposals that would be rejected are discarded
    without simulation; on the toy example this gives about 2.5 times more
    effective samples per second.}
  \item{collapse_mu}{if \code{TRUE} (default) the exchange step targets
    \eqn{p(\varphi \mid \Sigma, y)}, with the hierarchical mean integrated out:
    the prior ratio uses
    \eqn{\varphi \sim N(1_K \otimes \mu_0, I_K \otimes \Sigma + J_K \otimes
    \Sigma_0)}. \eqn{\mu} is then drawn from its full conditional given the new
    \eqn{\varphi} and \eqn{\Sigma}, and \eqn{\Sigma} given \eqn{\varphi} and
    \eqn{\mu}: a partially collapsed Gibbs sampler (van Dyk and Park, 2008)
    with the same stationary distribution. It removes the dependence between
    \eqn{\varphi} and \eqn{\mu} that slows the plain Gibbs scheme when the
    layers move together; on the toy example it gives about 1.7 times more
    effective samples for the same number of iterations.}
  \item{n_info}{networks simulated to estimate the Fisher information.}
  \item{n_newton}{moment-matching Newton steps from the penalised MPLE
    towards the posterior mode before the chains start:
    \eqn{\varphi \leftarrow \varphi + \{\mathrm{Cov}(s) + P\}^{-1}
    \{s_{obs} - E_\varphi(s) - P(\varphi - m)\}}.}
  \item{dispersion}{variance inflation of the starting points of chains 2,
    3, \dots}
  \item{seed}{random seed (each chain gets its own seed).}
  \item{verbose}{print progress.}
}
\details{
The initial proposal scale is \eqn{2.38/\sqrt{2d}} for a block of \eqn{d}
parameters: the usual random-walk scale, reduced by \eqn{\sqrt 2} because the
exchange step adds noise to the acceptance ratio with about the same variance
as the exact log-likelihood ratio.

Statistics enter in the rescaled form of \code{\link{ws_model}}; use
\code{coef(fit, scale = "ergm")} for coefficients on the \pkg{ergm} scale.
}
\value{
An object of class \code{"bergm_w"} (and \code{"bergm_ws"} for signed
networks): the model, the prior, the chains (draws
of \eqn{\varphi}, \eqn{\mu} and the lower triangle of \eqn{\Sigma},
acceptance indicators, final proposal scales and covariances), the settings and
the initialisation (starting value, initial proposal covariance and a fit
check of the starting value).
}
\references{
Christen, J. A. and Fox, C. (2005). Markov chain Monte Carlo using an
approximation. \emph{Journal of Computational and Graphical Statistics}, 14,
795--810.

Haario, H., Saksman, E. and Tamminen, J. (2001). An adaptive Metropolis
algorithm. \emph{Bernoulli}, 7, 223--242.

Murray, I., Ghahramani, Z. and MacKay, D. (2006). MCMC for doubly-intractable
distributions. \emph{Proceedings of UAI}.

van Dyk, D. A. and Park, T. (2008). Partially collapsed Gibbs samplers: theory
and methods. \emph{Journal of the American Statistical Association}, 103,
790--796.
}
\seealso{\code{\link{summary.bergm_w}}, \code{\link{plot.bergm_w}},
\code{\link{ws_predictive}}, \code{\link{ws_aux_check}}}
\examples{
data(toy_ws)
dat <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)
datu <- ws_data(toy_ws$weights, toy_ws$thresholds)
\donttest{
fs <- ~ L(~edges, ~(z & x_k)) + L(~edges, ~(!z & x_k)) +
  gwespL(0.5, fixed = TRUE, L.base = ~(z & x_k), Ls.path = list(~(z & x_k), ~(z & x_k)))
fit <- bergm_ws(dat, fs, scale = "gwesp", n_iter = 300, burnin = 100,
                n_chains = 2, seed = 1)
summary(fit)
fitu <- bergm_w(datu, ~ L(~edges, ~x_k) + L(~triangle, ~x_k),
                n_iter = 300, burnin = 100, n_chains = 2, seed = 1)
summary(fitu)
}
}
