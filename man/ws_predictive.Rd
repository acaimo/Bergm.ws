\name{ws_predictive}
\alias{ws_predictive}
\alias{print.ws_predictive}
\title{Prior and posterior predictive checks}
\description{
Simulates networks from the prior or from the posterior and compares them with
the observed network, both on the model statistics and on statistics that are
not in the model: for every layer, the numbers of triangles of each sign type
(\code{+++}, \code{-++}, \code{+--}, \code{---}) and the standard deviation of
the positive degrees.
}
\usage{
ws_predictive(object, type = c("posterior", "prior"), prior = NULL,
               nsim = 200, burnin = NULL, gof = TRUE, seed = NULL,
               verbose = TRUE)

\method{print}{ws_predictive}(x, digits = 2, ...)
}
\arguments{
  \item{object}{an \code{\link{bergm_w}} (posterior or prior check) or an
    \code{\link{ws_model}} (prior check only).}
  \item{type}{\code{"posterior"} or \code{"prior"}.}
  \item{prior}{prior for the prior check (default: that of the fit, or
    \code{ws_prior(model)}).}
  \item{nsim}{number of simulated networks, one per parameter draw.}
  \item{burnin}{MCMC steps per simulated network, started at the observed
    network; default five times the auxiliary chain length.}
  \item{gof}{also compute the statistics not in the model (needs the
    simulated networks, slower).}
  \item{seed}{random seed.}
  \item{verbose}{print progress.}
  \item{x}{an \code{ws_predictive} object.}
  \item{digits}{digits for printing.}
  \item{...}{ignored.}
}
\details{
The simulated networks are conditional on the observed lower layers (see
\code{\link{Bergm.ws-package}}), so each transition is checked given what is
observed below it. The mid-p value is the share of simulated values below
the observed one plus half the share of ties; values near 0 or 1 flag
statistics that the model does not reproduce.
}
\value{
An object of class \code{"ws_predictive"} with the simulated (\code{stats},
\code{gof}) and observed (\code{obs}, \code{gof_obs}) statistics and the
parameter draws used.
}
\seealso{\code{\link{bergm_w}}}
\examples{
data(toy_ws)

y <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)

f <- ~ L(~edges, ~(z & x_k)) + 
       L(~nodematch("pathway"), ~(z & x_k)) + L(~edges, ~(!z & x_k))

mod <- ws_model(y, f)
\donttest{
pp <- ws_predictive(mod, type = "prior", nsim = 50, seed = 1)
pp
}
}
