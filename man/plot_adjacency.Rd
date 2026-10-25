\name{plot_adjacency}
\alias{plot_adjacency}
\title{Signed adjacency matrix by weight level}
\description{
Heat map of the signed adjacency matrix: blue for positive and red for
negative edges, darker for edges that are still present in higher layers,
neutral grey for absent edges. Nodes are ordered by clustering the layout
coordinates or by an attribute, so that modules appear as blocks.
}
\usage{
plot_adjacency(data, order = c("layout", "attribute", "none"), by = NULL,
               labels = NULL, main = NULL)
}
\arguments{
  \item{data}{an \code{\link{ws_data}} object.}
  \item{order}{node order: hierarchical clustering of the
    \code{\link{ws_layout}} coordinates, a node attribute (\code{by}), or
    the original order.}
  \item{by}{attribute name for \code{order = "attribute"}.}
  \item{labels}{show node names (default when there are at most 50 nodes).}
  \item{main}{title.}
}
\value{The node order, invisibly.}
\examples{
data(toy_ws)
dat <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)
plot_adjacency(dat)
}
