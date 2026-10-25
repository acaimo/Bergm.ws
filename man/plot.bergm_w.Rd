\name{plot.bergm_w}
\alias{plot.bergm_w}
\title{Plots of a fitted model: traces, densities and posterior intervals}
\description{
Trace plots or posterior densities of the layer-specific parameters (one row
per layer, one column per statistic), the hierarchical means or the
between-layer standard deviations, with one colour per chain; or, with
\code{type = "intervals"}, posterior medians with 50\% and 95\% intervals of
every layer-specific parameter next to the hierarchical mean, one block per
statistic, so that the change of an effect across the layers and the amount of
pooling can be read at a glance.
}
\usage{
\method{plot}{bergm_w}(x, what = c("phi", "mu", "sd"),
     type = c("trace", "density", "intervals"), params = NULL,
     max_points = 2500, ...)
}
\arguments{
  \item{x}{an \code{bergm_w}.}
  \item{what}{which parameters (traces and densities).}
  \item{type}{\code{"trace"}, \code{"density"} or \code{"intervals"}.}
  \item{params}{optional subset of parameter names, as in
    \code{summary(x)}.}
  \item{max_points}{trace plots show at most this many evenly spaced draws
    per chain (all draws are used for densities and intervals).}
  \item{...}{for \code{type = "intervals"}: \code{main}, the title.}
}
\value{\code{x}, invisibly.}
\seealso{\code{\link{summary.bergm_w}}}
