\name{mixeddeepmm_grid}
\alias{mixeddeepmm_grid}
\title{Fit and compare a grid of mixed deep mixture models}
\description{Fit a supplied list of candidate architectures and compare observed log likelihood approximations and BIC.}
\usage{mixeddeepmm_grid(data, types = NULL, candidates, control = list())}
\arguments{
\item{data}{Data frame.}
\item{types}{Response type specification; see \code{mixeddeepmm}.}
\item{candidates}{List of named lists containing \code{family}, \code{heads} and \code{layers}.}
\item{control}{Control settings passed to each model.}
}
\value{List with comparison \code{table} and \code{fits}.}
\details{Monte Carlo likelihood error affects BIC uncertainty. BIC is an approximate fit criterion and does not itself measure clustering quality.}
