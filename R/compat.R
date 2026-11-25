# Namespace-local bridges required by the original DGMM SEM routines.
# Dependencies are requested at runtime so mixed models only need base R.
rmvnorm <- function(n, mean, sigma, ...) {
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("Continuous DGMM requires the optional 'mvtnorm' package.")
  mvtnorm::rmvnorm(n, mean, sigma, ...)
}
dmvnorm <- function(x, mean, sigma, log = FALSE, ...) {
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("Continuous DGMM requires the optional 'mvtnorm' package.")
  mvtnorm::dmvnorm(x, mean, sigma, log = log, ...)
}
`is.positive.definite` <- function(x, ...) {
  if (!requireNamespace("corpcor", quietly = TRUE))
    stop("Continuous DGMM requires the optional 'corpcor' package.")
  corpcor::is.positive.definite(x, ...)
}
`make.positive.definite` <- function(x, ...) {
  if (!requireNamespace("corpcor", quietly = TRUE))
    stop("Continuous DGMM requires the optional 'corpcor' package.")
  corpcor::make.positive.definite(x, ...)
}
