\name{mixeddeepmm_gower_silhouette}
\alias{mixeddeepmm_gower_silhouette}
\title{Mean Gower silhouette width for a cluster partition}
\usage{mixeddeepmm_gower_silhouette(clusters, distance)}
\arguments{
  \item{clusters}{One cluster label per observation; external true labels
  are not used.}
  \item{distance}{Precomputed \code{dist} or \code{dissimilarity} object,
  typically from \code{mixeddeepmm_gower_distance}.}
}
\description{Compute the average silhouette width via
\code{cluster::silhouette}. Returns \code{NA_real_} for partitions with
fewer than two or at least n nonempty groups. This score evaluates a
completed clustering, not an EM objective or a replacement for MCEM.}
\value{Numeric mean silhouette width, or \code{NA_real_}.}
\seealso{\code{mixeddeepmm_gower_distance}}
