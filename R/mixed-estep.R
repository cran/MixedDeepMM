.mdmm_paths <- function(k) {
  z <- as.matrix(expand.grid(lapply(k,seq_len), KEEP.OUT.ATTRS=FALSE))
  storage.mode(z) <- "integer"
  z
}

.mdmm_path_cache <- function(pars, s) {
  qx <- pars$qx; qt <- pars$qt; q <- qx+qt
  muH <- numeric(qx); alphaH <- numeric(qx)
  A <- matrix(0,qx,qt); psi <- numeric(qx)
  for (h in seq_along(pars$head)) {
    hh <- pars$head[[h]]; a <- s[h]; ix <- hh$idx
    muH[ix] <- hh$eta[[a]]
    alphaH[ix] <- hh$delta[[a]]
    A[ix,] <- hh$A[[a]]
    psi[ix] <- hh$psi[[a]]
  }
  b <- s[length(s)]
  muT <- pars$tail$eta[[b]]
  alphaT <- pars$tail$delta[[b]]
  sigmaT <- diag(pars$tail$psi[[b]], qt)
  mu <- c(muH+drop(A %*% muT), muT)
  alpha <- c(alphaH+drop(A %*% alphaT), alphaT)
  cov <- rbind(cbind(A %*% sigmaT %*% t(A)+diag(psi,qx), A %*% sigmaT),
               cbind(sigmaT %*% t(A), sigmaT))
  # M2DGMM observation-level MFA: y_C | x_C, S_C has its own
  # (intercept, loadings, diagonal residual covariance) for each S_C.
  # For single-head M1DGMM, there is one shared Gaussian link.
  measurement <- if (pars$two) s[1L] else 1L
  bb <- pars$measurement$b[[measurement]]
  LL <- pars$measurement$L[[measurement]]
  vv <- pars$measurement$psi[[measurement]]
  p <- length(bb)
  H <- matrix(0,p,q)
  if (p) H[,pars$xC] <- LL
  m <- bb+drop(H %*% mu)
  a <- drop(H %*% alpha)
  S <- if (p) H %*% cov %*% t(H)+diag(vv,p) else matrix(0,0,0)
  S_inv <- if (p) chol2inv(chol(.mdmm_spd(S))) else matrix(0,0,0)
  gain <- if (p) cov %*% t(H) %*% S_inv else matrix(0,q,0)
  V <- .mdmm_spd(cov - gain %*% H %*% cov)
  chV <- chol(V)
  w <- prod(mapply(function(pp,aa) pp[aa], pars$mixing, as.list(s)))
  list(mu=mu, alpha=alpha, sigma=S, inv=S_inv, meanC=m, alphaC=a,
       gain=gain, cholV=chV, logpi=log(pmax(w,1e-300)), q=q)
}

.mdmm_cont_posterior <- function(y, cache, family, nu) {
  p <- length(y)
  d <- if (p) y-cache$meanC else numeric(0)
  d2 <- if (p) drop(crossprod(d,cache$inv %*% d)) else 0
  chi <- if (p) drop(crossprod(cache$alphaC,cache$inv %*% cache$alphaC)) else 0
  if (family == "gaussian") {
    ld <- if (p) .mdmm_mvnorm_log(y,cache$meanC,cache$sigma) else 0
    return(list(logd=ld, shape=NA_real_, rate=NA_real_, chi=0, d2=d2))
  }
  sh <- (nu+p)/2
  rate <- (nu+d2)/2
  if (family == "t" || chi < 1e-12) {
    detlog <- if (p) determinant(.mdmm_spd(cache$sigma),logarithm=TRUE)$modulus[[1L]] else 0
    ld <- lgamma(sh)-lgamma(nu/2)-.5*(p*log(nu*pi)+detlog)-sh*log1p(d2/nu)
  } else {
    # P(y) under H~inverse-Gamma(nu/2,nu/2), y|H~N(m+H*a,H*S).
    # Parameterisation uses GIG(lambda=-sh, chi=nu+d2, psi=a'S^-1*a).
    x <- sqrt(chi*(nu+d2))
    lb <- log(pmax(besselK(x,nu=sh,expon.scaled=TRUE),.Machine$double.xmin))-x
    lognorm <- -0.5*(p*log(2*pi) + determinant(.mdmm_spd(cache$sigma),logarithm=TRUE)$modulus[[1L]])
    cross <- if (p) drop(crossprod(d,cache$inv %*% cache$alphaC)) else 0
    ld <- (nu/2)*log(nu/2)-lgamma(nu/2)+lognorm+cross +
      log(2) +(-sh/2)*(log(nu+d2)-log(chi))+lb
  }
  list(logd=as.numeric(ld),shape=sh,rate=rate,chi=chi,d2=d2)
}

.mdmm_draw_h <- function(M, posterior, family, nu) {
  if (family == "gaussian") return(rep(1,M))
  if (family == "t" || posterior$chi < 1e-12)
    return(1 / stats::rgamma(M,shape=posterior$shape,rate=posterior$rate))
  if (!requireNamespace("GIGrvg", quietly=TRUE))
    stop("Mixed skew-t fits require the optional GIGrvg package.")
  GIGrvg::rgig(M, lambda=-posterior$shape,chi=2*posterior$rate,psi=posterior$chi)
}

.mdmm_estep <- function(dat, pars, M, save_draws=TRUE) {
  paths <- .mdmm_paths(pars$k)
  n <- dat$n; K <- nrow(paths)
  logl <- matrix(-Inf,n,K)
  log_sd <- matrix(0,n,K)
  records <- if (save_draws) vector("list",n*K) else NULL
  caches <- lapply(seq_len(K),function(s) .mdmm_path_cache(pars, paths[s,]))
  for (s in seq_len(K)) {
    cc <- caches[[s]]
    for (i in seq_len(n)) {
      y <- if (dat$pC) dat$C[i,] else numeric(0)
      post <- .mdmm_cont_posterior(y,cc,pars$family,pars$nu)
      h <- .mdmm_draw_h(M,post,pars$family,pars$nu)
      h <- pmin(pmax(h,1e-8),1e8)
      z <- matrix(stats::rnorm(M*cc$q),M,cc$q) %*% cc$cholV
      dm <- if (dat$pC) y-cc$meanC else numeric(0)
      for (j in seq_len(M)) {
        z[j,] <- z[j,]*sqrt(h[j]) + cc$mu + h[j]*cc$alpha +
          drop(cc$gain %*% (dm-h[j]*cc$alphaC))
      }
      x <- z[,seq_len(pars$qx),drop=FALSE]
      ldis <- .mdmm_discrete_log(x,dat,i,pars)
      lw <- ldis-.mdmm_logsumexp(ldis)
      ww <- exp(lw)
      # Delta-method conditional MC variance of log(mean(exp(ldis))).
      vs <- exp(ldis-max(ldis))
      log_sd[i,s] <- if (M>1L) sqrt(stats::var(vs)/M)/pmax(mean(vs),1e-200) else NA_real_
      logl[i,s] <- cc$logpi+post$logd+.mdmm_logmeanexp(ldis)
      if (save_draws)
        records[[(s-1L)*n+i]] <- list(z=z,h=h,w=ww,i=i,s=s)
    }
  }
  norm <- apply(logl,1,.mdmm_logsumexp)
  tau <- exp(logl-norm)
  if (save_draws) for (s in seq_len(K)) for (i in seq_len(n)) {
    idx <- (s-1L)*n+i
    records[[idx]]$w <- records[[idx]]$w*tau[i,s]
  }
  # Estimated MC SE of total observed log likelihood under independent draws.
  loglik_se <- sqrt(sum((tau*log_sd)^2,na.rm=TRUE))
  list(loglik=sum(norm), mc_se=loglik_se,tau=tau,
       log_components=logl,paths=paths, records=records,
       first_layer=matrix(vapply(seq_len(pars$k[1]),function(a) rowSums(tau[,paths[,1]==a,drop=FALSE]),numeric(n)),nrow=n))
}
