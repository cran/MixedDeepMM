\name{simulate_mixeddeepmm}
\alias{simulate_mixeddeepmm}
\title{Simulate mixed observations and latent class labels}
\description{Simulate the manuscript's three-continuous and five-discrete response design, with three generating classes.}
\usage{simulate_mixeddeepmm(n = 600L, seed = 2026L)}
\arguments{\item{n}{Sample size (at least 30).}\item{seed}{Random seed.}}
\value{List with \code{data}, \code{types}, \code{group}, and \code{latent}.}
\examples{sim <- simulate_mixeddeepmm(n = 50, seed = 42)
str(sim$data)}
