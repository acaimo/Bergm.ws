\name{plot.ws_data}
\alias{plot.ws_data}
\alias{ws_layout}
\title{Plot a weighted network, signed or not}
\description{
Draws the baseline layer with positive edges in blue and negative edges in red
(dashed); without signs all edges are blue. The weight level of each edge is encoded by width and opacity, so
that the edges that survive to the highest layers stand out. Nodes are sized by
their weighted degree and can be coloured by an attribute (up to three
categories; further ones are shown in grey, and node shape gives a second
cue).

\code{ws_layout} computes the node positions. \code{type = "stress"} uses
stress majorisation on graph distances in which strong positive edges are short
and negative edges are long, so that cohesive positive groups sit together and
antagonistic nodes are pushed apart. \code{type = "circle"} places the nodes on
a circle, grouped by the attribute \code{by} (groups in the angular order of
the stress layout, nodes within a group in spectral order so that strongly tied
nodes are neighbours); edges are then bent towards the centre, every node is
labelled radially and each group is marked by an arc in its colour.
}
\usage{
\method{plot}{ws_data}(x, layout = NULL, color_by = NULL,
     labels = c("auto", "none", "all"), main = NULL, curvature = 0.12,
     neg_lty = "22", legend = TRUE, ...)

ws_layout(data, type = c("stress", "circle"), by = NULL,
           iterations = 500, neg_length = 1.5)
}
\arguments{
  \item{x, data}{an \code{\link{ws_data}} object.}
  \item{layout}{node coordinates from \code{ws_layout} (or any
    \eqn{n \times 2} matrix), or \code{"stress"} or \code{"circle"}; the
    default is the stress layout. With \code{"circle"}, nodes are grouped by
    \code{color_by}.}
  \item{color_by}{name of a node attribute for the node colours.}
  \item{labels}{\code{"auto"} labels every node in networks with up to 25
    nodes (80 in a circular layout) and the six with the largest weighted
    degree otherwise, each label placed where it overlaps the fewest nodes;
    or a logical vector.}
  \item{main}{title.}
  \item{curvature}{curvature of the edges (0 for straight lines).}
  \item{neg_lty}{line type of negative edges.}
  \item{legend}{draw the legend.}
  \item{type}{\code{"stress"} or \code{"circle"}.}
  \item{by}{for \code{type = "circle"}: name of the node attribute that
    defines the groups.}
  \item{iterations}{maximum number of stress-majorisation iterations.}
  \item{neg_length}{target length of negative edges, relative to 1 for the
    weakest positive edges.}
  \item{...}{ignored.}
}
\value{The layout, invisibly (\code{plot}); an \eqn{n \times 2} matrix with
attribute \code{"type"} (\code{ws_layout}).}
\seealso{\code{\link{plot_layers}}, \code{\link{plot_adjacency}}}
\examples{
data(toy_ws)
dat <- ws_data(toy_ws$weights, toy_ws$thresholds, toy_ws$signs, toy_ws$attributes)
xy <- ws_layout(dat)
plot(dat, layout = xy, color_by = "module")
plot(dat, layout = "circle", color_by = "module")
}
