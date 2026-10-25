\name{ws_aux_check}
\alias{ws_aux_check}
\alias{print.ws_aux_check}
\title{Are the auxiliary chains long enough?}
\description{
The exchange algorithm needs auxiliary networks drawn from the model. They are
obtained by running \code{aux_iters} MCMC steps from the observed network, so
they are only approximately exact. At the posterior mean, this function
compares the statistics of \code{n} such networks with those of \code{n}
networks from chains \code{factor} times longer.
}
\usage{
ws_aux_check(fit, n = 100, factor = 20, seed = NULL)

\method{print}{ws_aux_check}(x, ...)
}
\arguments{
  \item{fit}{an \code{\link{bergm_w}}.}
  \item{n}{networks per chain length.}
  \item{factor}{length of the long chains relative to \code{aux_iters}.}
  \item{seed}{random seed.}
  \item{x}{an \code{ws_aux_check} object.}
  \item{...}{ignored.}
}
\value{
A data frame with, for each statistic, the means of the short and long
chains, the long-chain standard deviation, the standardised difference of the
means and its z-score. A statistic is flagged when the difference exceeds a
quarter of a standard deviation and its z-score exceeds the Bonferroni
threshold \eqn{\Phi^{-1}(1 - 0.025/m)} for \eqn{m} statistics, so that a single
large z-score among many statistics is not mistaken for a problem; flagged
statistics suggest increasing \code{aux_iters}. With fewer than about 100
networks the Monte Carlo error of the means is large.
}
\seealso{\code{\link{bergm_w}}}
