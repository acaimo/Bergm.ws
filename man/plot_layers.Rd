\name{plot_layers}
\alias{plot_layers}
\title{Plot the nested layers side by side}
\description{
One panel per layer \eqn{x_0, \dots, x_K}, with a common layout, showing how
the network thins out and how the share of positive edges changes as the
weight threshold increases. Nodes without edges in a layer are greyed out.
}
\usage{
plot_layers(data, layout = NULL, color_by = NULL, ncol = NULL,
            curvature = 0.12, neg_lty = "22", main = NULL)
}
\arguments{
  \item{data}{an \code{\link{ws_data}} object.}
  \item{layout}{node coordinates from \code{\link{ws_layout}}, or
    \code{"stress"} (default) or \code{"circle"}.}
  \item{color_by}{name of a node attribute for the node colours.}
  \item{ncol}{number of panel columns (default: all layers in one row).}
  \item{curvature, neg_lty}{as in \code{\link{plot.ws_data}}.}
  \item{main}{overall title.}
}
\value{The layout, invisibly.}
\examples{
data(toy_ws)
dat <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)
plot_layers(dat, color_by = "module")
}
