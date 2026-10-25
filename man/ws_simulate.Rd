\name{ws_simulate}
\alias{ws_simulate}
\title{Simulate a small weighted network, signed or not}
\description{
Generates a synthetic weighted network with modular structure, for
examples and tests. Nodes belong to modules placed on a circle; the weight of
a dyad decreases with the distance between its nodes, so that strong edges are
clustered. With \code{signed = TRUE}, a minority faction has mostly negative and weak
edges with the other nodes, and signs are flipped more often for weak edges, so that the share
of positive edges increases with the weight.
}
\usage{
ws_simulate(n = 40, K = 3, n_modules = 3, density = 0.30,
            retention = c(0.55, 0.50, 0.50), signed = TRUE,
            faction_share = 0.18, sign_noise = 0.08, seed = NULL)
}
\arguments{
  \item{n}{number of nodes.}
  \item{K}{number of layer transitions (layers \eqn{x_0, \dots, x_K}).}
  \item{n_modules}{number of modules.}
  \item{density}{share of dyads in the baseline layer \eqn{x_0}.}
  \item{retention}{share of the edges of \eqn{x_{k-1}} kept in \eqn{x_k},
    recycled to length \code{K}.}
  \item{signed}{if \code{FALSE}, no signs and no faction effect.}
  \item{faction_share}{probability that a node belongs to the minority
    faction.}
  \item{sign_noise}{probability of flipping the sign of the weakest edges;
    it decreases linearly to zero for the strongest ones.}
  \item{seed}{optional random seed.}
}
\value{
A list with \code{weights} (symmetric matrix in (0, 1)), \code{signs}
(symmetric matrix of \eqn{\pm 1}, only if \code{signed}), \code{thresholds} (length \code{K + 1})
and \code{attributes} (data frame with \code{module} and \code{pathway}),
ready for \code{\link{ws_data}}.
}
\seealso{\code{\link{toy_ws}}, \code{\link{ws_data}}}
\examples{
sim <- ws_simulate(n = 30, seed = 1)
dat <- ws_data(sim$weights, sim$thresholds, sim$signs, sim$attributes)
dat
u <- ws_simulate(n = 30, signed = FALSE, seed = 1)
ws_data(u$weights, u$thresholds, attributes = u$attributes)
}
