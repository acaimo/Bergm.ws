\name{ws_model}
\alias{ws_model}
\alias{print.ws_model}
\title{Specify a multi-layer dissolution ERGM with an ergm.multi formula}
\description{
Takes an \pkg{ergm.multi} formula on the layers of a \code{\link{ws_data}}
object and checks it: every statistic must refer to exactly one layer, the
statistics must be grouped by layer, and every layer must have the same
statistics, because the hierarchical prior pools the \eqn{j}-th statistic of
all layers. Statistics equal to zero in the observed network and statistics
whose change statistics are nearly collinear on the dyads that can change are
reported (their parameters are then separated mostly by the prior).
}
\usage{
ws_model(data, formula, scale = NULL, labels = NULL)

\method{print}{ws_model}(x, ...)
}
\arguments{
  \item{data}{a \code{\link{ws_data}} object.}
  \item{formula}{a one-sided \pkg{ergm.multi} formula. The layers are called
    \code{x_1}, \dots, \code{x_K} and, for signed data, the signs \code{z}.
    A formula that uses \code{x_k} instead of numbered layers is a template:
    it is written once and repeated for \eqn{k = 1, \dots, K}.}
  \item{scale}{multipliers of the statistics: \code{NULL} (none),
    \code{"gwesp"} (every GWESP statistic with fixed decay \eqn{\alpha} is
    multiplied by \eqn{e^{-\alpha}}), or a numeric vector of length \eqn{p}
    (recycled over the layers) or \eqn{Kp}.}
  \item{labels}{optional names of the \eqn{p} statistics of a layer; by
    default they are derived from the \pkg{ergm} names, for example
    \code{edges+} for \code{L(~edges, ~(z & x_k))} and \code{gwesp-++} for a
    GWESP term on negative edges with positive shared partners.}
  \item{x}{a \code{ws_model} object.}
  \item{...}{ignored.}
}
\details{
Any \pkg{ergm} or \pkg{ergm.multi} term can be used, applied to a layer with
\code{L()} or with the layer-logic terms of \pkg{ergm.multi}. In signed data,
\code{(z & x_k)} selects the positive edges of layer \eqn{k} and
\code{(!z & x_k)} its negative edges; for example

\preformatted{
~ L(~edges, ~(z & x_k)) + L(~edges, ~(!z & x_k)) +
  gwespL(0.5, fixed = TRUE, L.base = ~(z & x_k),
         Ls.path = list(~(z & x_k), ~(z & x_k)))
}

models the persistence of positive and negative edges and positive triangle
closure in every layer. Without signs the layers are used directly, as in
\code{L(~edges, ~x_k)}. Objects named in the formula (a decay value, a vector
of attributes) are taken from its environment.

The coefficient of a statistic multiplied by \code{scale} is on the scale of
the multiplied statistic; with \code{scale = "gwesp"} a GWESP coefficient is
the change in log-odds for an edge with many shared partners, on the same
scale in every layer, which helps the hierarchical prior pool comparable
quantities. \code{coef(fit, scale = "ergm")} returns coefficients on the
original \pkg{ergm} scale.

A warning is given when the \eqn{j}-th statistics of different layers are
not the same term.
}
\value{
An object of class \code{"ws_model"} with the expanded formula (\code{rhs}),
the template if one was used, the statistic labels, the observed statistics
(\code{s_obs}, and \code{s_obs_scaled} after multiplication by \code{scale}),
the multipliers \code{scale} and the change-statistic correlations
\code{change_cor}.
}
\seealso{\code{\link{ws_prior}}, \code{\link{bergm_w}}}

