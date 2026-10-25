\name{ws_prior}
\alias{ws_prior}
\alias{print.ws_prior}
\title{Hierarchical prior for the layer-specific parameters}
\description{
Normal-inverse-Wishart prior:
\eqn{\varphi_k \mid \mu, \Sigma \sim N(\mu, \Sigma)} independently for
\eqn{k = 1, \dots, K}, \eqn{\mu \sim N(\mu_0, \Sigma_0)} and
\eqn{\Sigma \sim IW(\nu_0, S_0)}, with \eqn{E[\Sigma] = S_0 / (\nu_0 - p - 1)}.
}
\usage{
ws_prior(model, mu0 = 0, Sigma0 = 3, nu0 = NULL, S0 = 1)

\method{print}{ws_prior}(x, ...)
}
\arguments{
  \item{model}{an \code{\link{ws_model}} (or the number of statistics per
    layer).}
  \item{mu0}{prior mean of \eqn{\mu}, recycled.}
  \item{Sigma0}{prior covariance of \eqn{\mu}: a scalar (times the identity)
    or a matrix.}
  \item{nu0}{inverse-Wishart degrees of freedom; default \eqn{p + 9}, so that
    \eqn{E[\Sigma] = S_0 / 8}, a between-layer standard deviation of about
    0.35 when \eqn{S_0 = I}.}
  \item{S0}{inverse-Wishart scale: a scalar (times the identity) or a matrix.}
  \item{x}{an \code{ws_prior} object.}
  \item{...}{ignored.}
}
\details{
With few layers, the between-layer covariance \eqn{\Sigma} is determined
mostly by its prior: check the sensitivity of the conclusions to \code{nu0}
and \code{S0}.
}
\value{An object of class \code{"ws_prior"}.}
\examples{
data(toy_ws)
dat <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)
mod <- ws_model(dat, ~ L(~edges, ~(z & x_k)) + L(~edges, ~(!z & x_k)))
ws_prior(mod)
}
