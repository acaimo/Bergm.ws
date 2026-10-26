\name{ws_data}
\alias{ws_data}
\alias{print.ws_data}
\title{Nested layers from edge weights, with optional signs}
\description{
Thresholds a weighted network into nested binary layers
\eqn{x_k = \{ij : w_{ij} \ge t_k\}}, \eqn{k = 0, \dots, K}, optionally records the
sign of every edge of \eqn{x_0}, and builds the multilayer network and the
support constraints used by the model: \eqn{x_k} varies on the dyads of the
observed \eqn{x_{k-1}}, and the signs, when given, on the dyads of
\eqn{x_0}. Without signs the data are for \code{\link{bergm_w}}, with signs
for \code{\link{bergm_ws}}.
}
\usage{
ws_data(weights, thresholds, signs = NULL, attributes = NULL)

\method{print}{ws_data}(x, ...)
}
\arguments{
  \item{weights}{symmetric numeric matrix of edge weights (for example
    interaction confidence scores); \code{NA} means no edge.}
  \item{signs}{optional symmetric matrix whose sign gives the sign of each
    edge (for example a correlation matrix). Only dyads of \eqn{x_0} are used, and they
    must have a non-zero sign.}
  \item{thresholds}{strictly increasing vector of length \eqn{K + 1}; the
    first defines the baseline layer \eqn{x_0}.}
  \item{attributes}{optional data frame of node attributes, one row per node,
    available to the terms of the model formula.}
  \item{x}{an \code{ws_data} object.}
  \item{...}{ignored.}
}
\value{
An object of class \code{"ws_data"}: a list with the multilayer network
\code{y} (signs, if any, and layers \eqn{x_1, \dots, x_K}), the layers \code{x}
(\eqn{x_0, \dots, x_K}), the matrix \code{levels} (for each dyad, the highest
layer containing it, or -1), the sign matrix \code{signs}, the support
\code{constraint} formula and the inputs.
}
\seealso{\code{\link{ws_model}}}
\examples{
data(toy_ws)

y <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)
y

x <- ws_data(toy_ws$weights, toy_ws$thresholds)
x
}
