.mdmm_logsumexp <- function(x) {
  z <- max(x)
  if (!is.finite(z)) return(-Inf)
  z + log(sum(exp(x - z)))
}
.mdmm_softmax <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(rep(1 / length(x), length(x)))
  y <- exp(x - m)
  y / sum(y)
}
.mdmm_spd <- function(a, floor = 1e-7) {
  a <- (a + t(a)) / 2
  e <- eigen(a, symmetric = TRUE)
  tcrossprod(e$vectors %*% diag(sqrt(pmax(e$values, floor)), nrow = length(e$values)))
}
.mdmm_mvnorm_log <- function(y, mean, sigma) {
  if (!length(mean)) return(0)
  ch <- chol(.mdmm_spd(sigma))
  v <- backsolve(ch, y - mean, transpose = TRUE)
  -0.5 * (length(mean) * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(v^2))
}
.mdmm_weighted_ls <- function(x, y, weights, ridge = 1e-6) {
  # y matrix, x matrix, with ridge to handle sparse mixture paths.
  sw <- sqrt(pmax(weights, 0))
  xw <- x * sw
  yw <- y * sw
  p <- ncol(x)
  reg <- diag(ridge, p)
  if (p) reg[1, 1] <- 0
  solve(crossprod(xw) + reg + diag(1e-9, p), crossprod(xw, yw))
}
.mdmm_prob <- function(eta) plogis(pmax(-35, pmin(35, eta)))
.mdmm_logmeanexp <- function(logx) .mdmm_logsumexp(logx) - log(length(logx))
.mdmm_finite <- function(x, name) {
  if (!is.numeric(x) || any(!is.finite(x))) stop(name, " must be finite numeric data.")
}
