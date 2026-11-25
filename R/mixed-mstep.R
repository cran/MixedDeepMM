.mdmm_stack_draws <- function(estep, pars) {
  rr <- estep$records
  nrec <- length(rr)
  N <- sum(vapply(rr,function(z)length(z$w),integer(1)))
  z <- matrix(0,N,pars$qx+pars$qt)
  h <- w <- numeric(N)
  obs <- path <- integer(N)
  pos <- 1L
  for (j in seq_len(nrec)) {
    a <- rr[[j]]; m <- length(a$w)
    ix <- seq.int(pos,length.out=m)
    z[ix,] <- a$z; h[ix] <- a$h; w[ix] <- a$w
    obs[ix] <- a$i; path[ix] <- a$s
    pos <- pos+m
  }
  list(z=z,h=h,w=w,obs=obs,path=path)
}

.mdmm_opt_binary_count <- function(X, y, weights, current, type, ridge) {
  f <- function(b) {
    eta <- drop(X %*% b)
    if (type == "binary") {
      ll <- sum(weights*(y*eta - pmax(eta,0)-log1p(exp(-abs(eta)))))
    } else {
      ll <- sum(weights*(y*eta-exp(pmin(eta,20))-lgamma(y+1)))
    }
    -ll + .5*ridge*sum(b[-1L]^2)
  }
  g <- function(b) {
    eta <- drop(X %*% b)
    m <- if (type == "binary") .mdmm_prob(eta) else exp(pmin(eta,20))
    grad <- -drop(crossprod(X,weights*(y-m)))
    grad[-1L] <- grad[-1L]+ridge*b[-1L]
    grad
  }
  fit <- try(stats::optim(current,f,g,method="L-BFGS-B",lower=-12,upper=12,
                          control=list(maxit=100)),silent=TRUE)
  if (inherits(fit,"try-error") || !all(is.finite(fit$par))) return(current)
  fit$par
}

.mdmm_opt_nominal <- function(X, y, weights, current, ridge) {
  levels <- nrow(current); nx <- ncol(X)
  objective <- function(b, gradient=FALSE) {
    B <- matrix(b,levels,nx)
    eta <- cbind(X %*% t(B),0)
    mx <- apply(eta,1,max)
    ee <- exp(eta-mx)
    pp <- ee/rowSums(ee)
    if (!gradient) {
      ll <- sum(weights*log(pmax(pp[cbind(seq_len(nrow(X)),y)],1e-300)))
      return(-ll + .5*ridge*sum(B[,-1L,drop=FALSE]^2))
    }
    resid <- -pp[,seq_len(levels),drop=FALSE]
    for (j in seq_len(levels)) resid[,j] <- resid[,j]+as.numeric(y==j)
    grad <- -t(resid*weights) %*% X
    grad[,-1L] <- grad[,-1L]+ridge*B[,-1L]
    as.vector(grad)
  }
  fit <- try(stats::optim(as.vector(current), objective,
                  gr=function(b) objective(b,TRUE),method="L-BFGS-B",
                  lower=-12,upper=12,control=list(maxit=100)),silent=TRUE)
  if (inherits(fit,"try-error") || !all(is.finite(fit$par))) return(current)
  matrix(fit$par,levels,nx)
}

.mdmm_opt_ordinal <- function(X, y, weights, current, ridge) {
  m <- length(current$thresholds)
  # theta[1] free; differences strictly positive, parameterised on log scale.
  old <- c(current$thresholds[1],log(pmax(diff(current$thresholds),1e-3)),current$beta)
  objective <- function(a) {
    theta <- if (m > 1L) c(a[1],a[1]+cumsum(exp(a[2:m]))) else a[1]
    beta <- tail(a,ncol(X))
    eta <- drop(X %*% beta)
    cp <- sapply(theta,function(th) .mdmm_prob(th-eta))
    cp <- matrix(cp,nrow=nrow(X))
    p <- cbind(cp[,1], if(m>1L) cp[,2:m,drop=FALSE]-cp[,1:(m-1L),drop=FALSE],1-cp[,m])
    -sum(weights*log(pmax(p[cbind(seq_len(nrow(X)),y)],1e-200))) + .5*ridge*sum(beta^2)
  }
  lower <- c(-12, rep(-6,max(0L,m-1L)),rep(-12,ncol(X)))
  upper <- c(12, rep(6,max(0L,m-1L)),rep(12,ncol(X)))
  fit <- try(stats::optim(old,objective,method="L-BFGS-B",lower=lower,upper=upper,
                         control=list(maxit=65)),silent=TRUE)
  if (inherits(fit,"try-error") || !all(is.finite(fit$par))) return(current)
  a <- fit$par
  list(type="ordinal", thresholds=if (m > 1L) c(a[1],a[1]+cumsum(exp(a[2:m]))) else a[1],
       beta=tail(a,ncol(X)))
}

.mdmm_update <- function(dat, pars, estep, ridge=1e-3, floor=.03, estimate_nu=FALSE) {
  rec <- .mdmm_stack_draws(estep,pars)
  paths <- estep$paths
  N <- nrow(rec$z)
  tlat <- rec$z[,pars$qx+seq_len(pars$qt),drop=FALSE]
  h <- rec$h; u <- 1/h
  w <- rec$w
  uu <- w*u
  out <- pars
  # Independent allocation priors: count effective marginal responsibilities.
  for (j in seq_along(pars$k)) {
    pp <- vapply(seq_len(pars$k[j]),function(a)sum(estep$tau[,paths[,j]==a,drop=FALSE]),numeric(1))
    out$mixing[[j]] <- (pp+.01)/sum(pp+.01)
  }
  # Weighted Gaussian regressions for x_head | t,H,s_head.
  for (j in seq_along(out$head)) {
    hh <- out$head[[j]]
    for (a in seq_len(pars$k[j])) {
      ix <- which(paths[rec$path,j] == a)
      if (length(ix) < 2L || sum(w[ix]) < 2) next
      X <- cbind(1,tlat[ix,,drop=FALSE])
      if (pars$family == "skew-t") X <- cbind(X,h[ix])
      Y <- rec$z[ix,hh$idx,drop=FALSE]
      B <- try(.mdmm_weighted_ls(X,Y,uu[ix],ridge),silent=TRUE)
      if (inherits(B,"try-error") || any(!is.finite(B))) next
      hh$eta[[a]] <- as.numeric(B[1,])
      hh$A[[a]] <- t(B[1+seq_len(pars$qt),,drop=FALSE])
      if (pars$family == "skew-t") hh$delta[[a]] <- as.numeric(B[nrow(B),])
      res <- Y-X %*% B
      hh$psi[[a]] <- pmax(floor,colSums(res^2*uu[ix])/sum(w[ix]))
    }
    out$head[[j]] <- hh
  }
  # Gaussian-mixture shared-tail latent distribution, equivalently the
  # marginal of a final Gaussian factor layer for a one-dimensional tail.
  # The deeper Gaussian root is integrated out, not estimated separately.
  j <- length(pars$k)
  for (a in seq_len(pars$k[j])) {
    ix <- which(paths[rec$path,j] == a)
    if (length(ix)<2L || sum(w[ix]) < 2) next
    mu <- colSums(tlat[ix,,drop=FALSE]*uu[ix])/sum(uu[ix])
    out$tail$eta[[a]] <- as.numeric(mu)
    residual <- sweep(tlat[ix,,drop=FALSE],2L,mu,'-')
    out$tail$psi[[a]] <- pmax(floor,colSums(residual^2*uu[ix])/sum(w[ix]))
  }
  # M2DGMM: component-specific MFA observation regressions; M1DGMM:
  # a common Gaussian GLLVM link (diagonal conditional residual covariance).
  if (dat$pC) {
    nc <- length(out$measurement$b)
    for (a in seq_len(nc)) {
      ix <- if (pars$two) which(paths[rec$path,1L]==a) else seq_len(N)
      if (length(ix)<2L || sum(w[ix])<2) next
      X <- cbind(1,rec$z[ix,pars$xC,drop=FALSE])
      Y <- dat$C[rec$obs[ix],,drop=FALSE]
      B <- try(.mdmm_weighted_ls(X,Y,uu[ix],ridge),silent=TRUE)
      if (inherits(B,'try-error') || any(!is.finite(B))) next
      out$measurement$b[[a]] <- as.numeric(B[1,])
      out$measurement$L[[a]] <- t(B[-1L,,drop=FALSE])
      res <- Y-X %*% B
      out$measurement$psi[[a]] <- pmax(floor,colSums(res^2*uu[ix])/sum(w[ix]))
    }
  }
  # Discrete conditional GLMs: weights are tau * conditional IS weights,
  # NOT precision-weighted likelihood contributions.
  XD <- rec$z[,pars$xD,drop=FALSE]
  X <- cbind(1,XD)
  for (nm in names(pars$links)) {
    link <- pars$links[[nm]]
    y <- dat$D[[nm]][rec$obs]
    if (link$type %in% c("binary","count")) {
      link$beta <- .mdmm_opt_binary_count(X,y,w,link$beta,link$type,ridge)
    } else if (link$type == "nominal") {
      link$beta <- .mdmm_opt_nominal(X,y,w,link$beta,ridge)
    } else if (link$type == "ordinal") {
      link <- .mdmm_opt_ordinal(XD,y,w,link,ridge)
    }
    out$links[[nm]] <- link
  }
  if (estimate_nu && pars$family != "gaussian") {
    mean_term <- sum(w*(-log(h)-1/h))/sum(w)
    root <- function(nu) abs(log(nu/2)-digamma(nu/2)+1+mean_term)
    out$nu <- stats::optimize(root,interval=c(2.1,120))$minimum
  }
  .mdmm_identify(out)
}

.mdmm_parameters <- function(dat,pars,estimate_nu=FALSE) {
  k <- pars$k; qt <- pars$qt; qd <- length(pars$xD)
  np <- sum(k-1L)
  for (h in seq_along(pars$head)) {
    q <- pars$head[[h]]$q
    np <- np+k[h]*(q+q*qt+q+if(pars$family=="skew-t") q else 0L)
  }
  np <- np+tail(k,1L)*(2L*qt) # tail means and diagonal covariances
  np <- np+length(pars$measurement$b)*dat$pC*(2L+length(pars$xC))
  for (nm in names(pars$links)) {
    li <- pars$links[[nm]]
    if (li$type %in% c("binary","count")) np <- np+qd+1L
    if (li$type=="ordinal") np <- np+length(li$thresholds)+qd
    if (li$type=="nominal") np <- np+nrow(li$beta)*(qd+1L)
  }
  if (estimate_nu && pars$family!="gaussian") np <- np+1L
  # First tail/head locations are reference-level constraints.
  np <- np-qt-sum(vapply(pars$head,`[[`,integer(1),"q"))
  np
}

# Anchor the first tail and head location to zero. These transformations leave
# the implied probability model unchanged, unlike hard zeroing which would not.
.mdmm_identify <- function(pars) {
  ref <- pars$tail$eta[[1L]]
  for (a in seq_along(pars$tail$eta))
    pars$tail$eta[[a]] <- pars$tail$eta[[a]] - ref
  for (h in seq_along(pars$head)) {
    head <- pars$head[[h]]
    for (a in seq_along(head$eta))
      head$eta[[a]] <- head$eta[[a]] + drop(head$A[[a]] %*% ref)
    pars$head[[h]] <- head
  }
  for (h in seq_along(pars$head)) {
    head <- pars$head[[h]]
    ref <- head$eta[[1L]]
    for (a in seq_along(head$eta)) head$eta[[a]] <- head$eta[[a]] - ref
    pars$head[[h]] <- head
    if (!pars$two || h==1L) {
      for (a in seq_along(pars$measurement$b))
        if (length(pars$measurement$b[[a]]))
          pars$measurement$b[[a]] <- pars$measurement$b[[a]] +
            drop(pars$measurement$L[[a]] %*% ref)
    }
    if (!pars$two || h==2L) {
      for (nm in names(pars$links)) {
        li <- pars$links[[nm]]
        if (li$type == "ordinal") {
          li$thresholds <- li$thresholds - sum(li$beta*ref)
        } else if (li$type == "nominal") {
          li$beta[,1] <- li$beta[,1] + drop(li$beta[,-1L,drop=FALSE] %*% ref)
        } else {
          li$beta[1] <- li$beta[1] + sum(li$beta[-1L]*ref)
        }
        pars$links[[nm]] <- li
      }
    }
  }
  pars
}
