\name{predict.mixeddeepmm}
\alias{predict.mixeddeepmm}
\alias{fitted.mixeddeepmm}
\title{Predict pathway posteriors, clusters and responses}
\description{Obtain posterior probabilities of latent pathways, their marginal layer probabilities, classifications or model-implied response means/probabilities.}
\usage{
\method{predict}{mixeddeepmm}(
  object,
  newdata = NULL,
  type = c("cluster", "posterior", "path", "layer", "response"),
  layer = NULL,
  mc_samples = NULL,
  cluster_level = c(
    "paper", "first", "shared_tail", "continuous_head", "discrete_head"
  ),
  ...
)
\method{fitted}{mixeddeepmm}(object, ...)
}
\arguments{
\item{object}{Fitted \code{mixeddeepmm} object.}
\item{newdata}{Optional data frame with the same named columns and known categorical levels.}
\item{type}{Return specified-layer MAP classification (\code{cluster}), complete-path posterior (\code{posterior}), MAP full pathway (\code{path}), layer marginal posterior (\code{layer}), or marginal link probabilities and expected counts (\code{response}).}
\item{layer}{Optional integer layer index for \code{type="layer"} or \code{type="cluster"}; overrides cluster_level.}
\item{mc_samples}{Override MC importance sample size when predicting mixed data.}
\item{cluster_level}{Named layer when layer is NULL: paper uses the common first layer for single head and shared tail for two heads; first selects layer 1; shared_tail selects the last layer. Head options are available as indicated by the architecture.}
\item{...}{Additional arguments.}
}
\details{By default cluster_level=\code{"paper"} selects the first common head in single-head models and the shared tail in two-head models. When predicting link probabilities, observed responses inform the latent posterior. Thus the estimates are conditional fitted probabilities and are not leave-one-out predictive scores.}
\value{Numeric vector, matrix or named list of fitted values, according to \code{type}.}
