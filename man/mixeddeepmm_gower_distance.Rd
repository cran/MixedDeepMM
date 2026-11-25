\name{mixeddeepmm_gower_distance}
\alias{mixeddeepmm_gower_distance}
\title{Mixed-response Gower dissimilarity}
\usage{mixeddeepmm_gower_distance(data, types)}
\arguments{
  \item{data}{Complete mixed-response predictors, without the external truth.}
  \item{types}{Named character vector defining continuous, count, binary,
  nominal and ordinal response types.}
}
\description{Construct a single mixed-type dissimilarity matrix with equal
variable weights using \code{cluster::daisy(metric="gower")}. Binary and
nominal predictors are symmetric categorical variables; ordinal predictors
are ordered; continuous/count predictors are range-normalised numeric.}
\value{A \code{dist}/\code{dissimilarity} object.}
\seealso{\code{mixeddeepmm_gower_silhouette}}
