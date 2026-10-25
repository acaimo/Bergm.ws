\name{summary.bergm_w}
\alias{summary.bergm_w}
\alias{print.summary.bergm_w}
\alias{print.bergm_w}
\alias{coef.bergm_w}
\alias{as.mcmc.list.bergm_w}
\title{Posterior summaries and convergence diagnostics}
\description{
Posterior means, standard deviations and quantiles of the hierarchical means
\eqn{\mu}, the layer-specific parameters \eqn{\varphi_k} and the between-layer
standard deviations \eqn{\Sigma_{rr}^{1/2}}, with the potential scale
reduction factor \eqn{\hat R} and the effective sample size of every
parameter, the acceptance rate of each chain after burn-in and the largest
posterior correlation between parameters of different layers.
}
\usage{
\method{summary}{bergm_w}(object, probs = c(0.025, 0.5, 0.975), ...)

\method{print}{summary.bergm_w}(x, digits = 3, ...)

\method{print}{bergm_w}(x, ...)

\method{coef}{bergm_w}(object, what = c("phi", "mu"), scale = c("model", "ergm"), ...)

\method{as.mcmc.list}{bergm_w}(x, what = c("phi", "mu", "sd"), ...)
}
\arguments{
  \item{object, x}{an \code{bergm_w} (or its summary).}
  \item{probs}{posterior quantiles to report.}
  \item{digits}{digits for printing.}
  \item{what}{\code{"phi"}: layer-specific parameters; \code{"mu"}:
    hierarchical means; \code{"sd"}: between-layer standard deviations.}
  \item{scale}{\code{"model"}: the (possibly rescaled) scale of the model;
    \code{"ergm"}: coefficients of the unscaled \pkg{ergm} statistics.}
  \item{...}{ignored.}
}
\details{
\eqn{\hat R} uses \code{coda::gelman.diag} without discarding further draws
and needs at least two chains; values above 1.05, or effective sample sizes
below 100, call for longer chains. A large correlation between parameters of
different layers suggests \code{update = "joint"}.
}
\value{
\code{summary} returns an object of class \code{"summary.bergm_w"};
\code{coef} a \eqn{K \times p} matrix (or a vector for \code{what = "mu"});
\code{as.mcmc.list} a \code{coda} \code{mcmc.list}.
}
\seealso{\code{\link{bergm_w}}, \code{\link{plot.bergm_w}}}
