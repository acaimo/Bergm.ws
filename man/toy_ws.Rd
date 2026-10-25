\name{toy_ws}
\alias{toy_ws}
\docType{data}
\title{A small simulated signed weighted network}
\description{
A signed weighted network on 40 nodes simulated with
\code{ws_simulate(n = 40, seed = 33)}. It is entirely synthetic: three
modules of nodes with strong positive edges inside them, a minority faction
whose edges with the others tend to be negative and weak, and signs that are
noisier for weak edges. As in co-expression networks built from interaction
confidence scores, the share of positive edges grows with the weight (77\%,
87\%, 94\% and 97\% in the four layers).
}
\usage{data(toy_ws)}
\format{
A list with components
\describe{
  \item{\code{weights}}{40 x 40 symmetric matrix of edge weights in (0, 1).}
  \item{\code{signs}}{40 x 40 symmetric matrix of signs (\eqn{\pm 1}).}
  \item{\code{thresholds}}{the four thresholds defining layers \eqn{x_0}
    to \eqn{x_3}.}
  \item{\code{attributes}}{data frame with the node attributes \code{module}
    (three modules) and \code{pathway} (a binary annotation that does not
    affect the edges).}
}
}
\examples{
data(toy_ws)
dat <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)
dat
ws_data(toy_ws$weights, toy_ws$thresholds, attributes = toy_ws$attributes)
}
\keyword{datasets}
