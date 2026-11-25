\name{summary.mixeddeepmm}
\alias{summary.mixeddeepmm}
\alias{print.mixeddeepmm}
\alias{print.summary.mixeddeepmm}
\alias{logLik.mixeddeepmm}
\alias{plot.mixeddeepmm}
\title{Summary, print, likelihood and plot methods for deep mixtures}
\description{Display goodness of fit, Monte Carlo precision, iteration count and likelihood history (without an end-of-run convergence flag).}
\usage{
\method{summary}{mixeddeepmm}(object, ...)
\method{print}{mixeddeepmm}(x, ...)
\method{print}{summary.mixeddeepmm}(x, ...)
\method{logLik}{mixeddeepmm}(object, ...)
\method{plot}{mixeddeepmm}(x, ...)
}
\arguments{
\item{object}{A \code{mixeddeepmm} fit.}
\item{x}{A model or summary object.}
\item{...}{Other options passed to plot, as appropriate.}
}
\value{A summary list, an invisible object for print/plot, or \code{logLik} object.}
