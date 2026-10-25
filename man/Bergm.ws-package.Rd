\name{Bergm.ws-package}
\alias{Bergm.ws-package}
\alias{Bergm.ws}
\docType{package}
\title{Bayesian multi-layer dissolution ERGMs for weighted and signed weighted networks}
\description{
A weighted network is decomposed into nested binary layers
\eqn{x_0 \supseteq x_1 \supseteq \dots \supseteq x_K} by thresholding its edge
weights. Each transition \eqn{x_{k-1} \to x_k} is modelled by an exponential
random graph model whose statistics are computed on layer \eqn{k},
conditionally on the observed lower layers: \code{\link{bergm_w}} fits the
model to the layers alone. When the sign \eqn{z} of every edge of \eqn{x_0}
is also recorded, \code{\link{bergm_ws}} models the signs jointly with the
layers:
\deqn{p(z, x_{1:K} \mid \tilde x_{0:K-1}, \varphi_{1:K}) \propto
\exp\Big\{\sum_{k=1}^K \varphi_k^\top s(z, x_k)\Big\}
\prod_{k=1}^K 1\{x_k \subseteq \tilde x_{k-1}\},}{
p(z, x_1:K | observed x_0:K-1, phi) proportional to
exp{ sum_k phi_k' s(z, x_k) } prod_k 1{x_k subset of observed x_(k-1)},}
where \eqn{\tilde x_{k-1}} is the \emph{observed} preceding layer. The model is
therefore conditional on the observed lower layers. The layer-specific
parameters are partially pooled by a Normal-inverse-Wishart prior,
\eqn{\varphi_k \mid \mu, \Sigma \sim N(\mu, \Sigma)},
\eqn{\mu \sim N(\mu_0, \Sigma_0)}, \eqn{\Sigma \sim IW(\nu_0, S_0)}, and the
posterior is sampled with an adaptive approximate exchange algorithm
(\code{\link{bergm_w}}, \code{\link{bergm_ws}}): delayed acceptance with a Gaussian approximation of
the likelihood, so that most proposals that would be rejected are discarded
without simulating a network; the hierarchical mean integrated out of the
exchange step; and a proposal covariance initialised from the Fisher
information, learned during the burn-in, pooled across chains and frozen for
the retained draws.
}
\details{
The typical workflow is:
\enumerate{
  \item \code{\link{ws_data}}: build the layers from a weight matrix, a sign
    matrix and the thresholds (\code{\link{ws_simulate}} creates a toy
    example; \code{\link{toy_ws}} is one).
  \item \code{\link{plot.ws_data}}, \code{\link{plot_layers}},
    \code{\link{plot_adjacency}}: look at the network.
  \item \code{\link{ws_model}}: write the model as an \pkg{ergm.multi}
    formula on the layers \code{x_1}, \dots, \code{x_K} and the signs
    \code{z}, or once as a template with \code{x_k}. The formula is checked:
    every statistic must use exactly one layer and every layer the same
    statistics; nearly collinear statistics are flagged.
  \item \code{\link{ws_prior}}, \code{\link{ws_predictive}} with
    \code{type = "prior"}: prior and prior predictive check.
  \item \code{\link{bergm_w}} or \code{\link{bergm_ws}}: adaptive exchange sampler, then
    \code{\link{summary.bergm_w}} and \code{\link{plot.bergm_w}} for estimates
    and convergence diagnostics.
  \item \code{\link{ws_predictive}}: posterior predictive check, including
    statistics not used in the model; \code{\link{ws_aux_check}}: are the
    auxiliary chains long enough?
}
}
\seealso{
\code{vignette("Bergm.ws", package = "Bergm.ws")}
}
\references{
Caimo, A. and Gollini, I. (2020). A multi-layer exponential random graph
modelling approach for weighted networks. \emph{Computational Statistics &
Data Analysis}, 142, 106825.

Caimo, A. and Friel, N. (2011). Bayesian inference for exponential random
graph models. \emph{Social Networks}, 33, 41--55.

Murray, I., Ghahramani, Z. and MacKay, D. (2006). MCMC for doubly-intractable
distributions. \emph{Proceedings of UAI}.

Christen, J. A. and Fox, C. (2005). Markov chain Monte Carlo using an
approximation. \emph{Journal of Computational and Graphical Statistics}, 14,
795--810.

Haario, H., Saksman, E. and Tamminen, J. (2001). An adaptive Metropolis
algorithm. \emph{Bernoulli}, 7, 223--242.

Krivitsky, P. N., Koehly, L. M. and Marcum, C. S. (2020). Exponential-family
random graph models for multi-layer networks. \emph{Psychometrika}, 85,
630--659.

van Dyk, D. A. and Park, T. (2008). Partially collapsed Gibbs samplers: theory
and methods. \emph{Journal of the American Statistical Association}, 103,
790--796.
}
\author{Alberto Caimo, Paolo Dalena, Isabella Gollini}
\keyword{package}
