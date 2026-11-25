# Deterministic, data-derived initialisation for mixed-response MCEM.
# PCA/indicator embeddings, deterministic Lloyd clustering, Gaussian least
# squares and regularised response regressions; no random draws are made.
# This is NSEP-inspired but not an implementation of the published NSEP.
# Training labels must not appear in dat, as for any other model fitting step.

.mdmm_init_standardize <- function(x) {
  x <- as.matrix(x)
  if (!ncol(x)) return(matrix(0, nrow(x), 0L))
  center <- colMeans(x)
  s <- apply(x, 2L, stats::sd)
  s[!is.finite(s) | s < 1e-8] <- 1
  a <- sweep(sweep(x, 2L, center, '-'), 2L, s, '/')
  pmax(pmin(a, 5), -5) # matrix argument first: retain dimensions
}

.mdmm_init_discrete_embedding <- function(dat) {
  if (!dat$pD) return(matrix(0, dat$n, 0L))
  blocks <- lapply(names(dat$D), function(nm) {
    v <- dat$D[[nm]]
    type <- dat$types[[nm]]
    if (type == 'nominal') {
      nlev <- length(dat$schema[[nm]]$levels)
      indicators <- sapply(seq_len(nlev), function(j) as.numeric(v == j))
      indicators <- matrix(indicators, nrow=dat$n)
      return(.mdmm_init_standardize(indicators) / sqrt(nlev))
    }
    if (type == 'count') v <- log1p(v)
    .mdmm_init_standardize(matrix(as.numeric(v), ncol=1L))
  })
  do.call(cbind, blocks)
}

.mdmm_init_scores <- function(x, q) {
  x <- as.matrix(x)
  n <- nrow(x)
  z <- matrix(0, n, q)
  if (!ncol(x)) return(z)
  x <- sweep(x, 2L, colMeans(x), '-')
  # svd() is deterministic for fixed inputs, apart from irrelevant factor
  # sign conventions; resolve the sign so that equal inputs yield equal scores.
  sv <- svd(x, nu=min(n, q, ncol(x)), nv=0L)
  rank <- sum(is.finite(sv$d) & sv$d > max(sv$d,0)*1e-9)
  keep <- min(q, rank)
  if (keep) {
    z[,seq_len(keep)] <- sweep(sv$u[,seq_len(keep),drop=FALSE],
                               2L,sv$d[seq_len(keep)],'*')
    for (j in seq_len(keep)) {
      ix <- which.max(abs(z[,j]))
      if (z[ix,j] < 0) z[,j] <- -z[,j]
    }
    z[,seq_len(keep)] <- .mdmm_init_standardize(z[,seq_len(keep),drop=FALSE])
  }
  # Insufficient observed rank: unused latent columns stay at zero. They are
  # protected by positive covariance floors, not fabricated random signal.
  z
}

.mdmm_init_groups <- function(x, k, variant=1L) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (k > n) stop('Cannot initialise more components than observations.')
  if (k == 1L) return(rep(1L, n))
  q <- ncol(x)
  if (!q) stop('Cannot cluster a latent space with no columns.')
  # Different starts are deterministic alternative projections of the data.
  # For a one-dimensional tail, quantile offset also yields distinct splits.
  v <- as.numeric(variant)
  j <- seq_len(q)
  weights <- cos((v %% 37 + 1) * j * .173) + sin((v %% 23 + 1)*j*.271)
  if (sum(abs(weights)) < 1e-8) weights[1L] <- 1
  projection <- drop(x %*% weights)
  ord <- order(projection, seq_len(n))
  offset <- ((as.integer(v %% 9) - 4L) / 40) # [-.1, .1], no RNG
  probs <- pmax(.01,pmin(.99,((seq_len(k)-.5)/k)+offset))
  centers_idx <- ord[pmax(1L,pmin(n,as.integer(round(1+(n-1)*probs))))]
  center_mat <- x[centers_idx,,drop=FALSE]
  unique_centers <- !duplicated(as.data.frame(center_mat))
  if (sum(unique_centers) == k) {
    fit <- try(stats::kmeans(x, centers=center_mat, iter.max=50L,
                             algorithm='Lloyd'), silent=TRUE)
    if (!inherits(fit, 'try-error') && length(unique(fit$cluster))==k)
      return(as.integer(fit$cluster))
  }
  # Entirely deterministic, balanced partition for tied/rank-deficient data.
  grp <- integer(n)
  grp[ord] <- as.integer(cut(seq_len(n), breaks=k, labels=FALSE))
  grp
}

.mdmm_init_ridge <- function(X, Y, penalty=.05) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  pen <- diag(penalty,ncol(X)); pen[1L,1L] <- 1e-8
  B <- try(solve(crossprod(X)+pen,crossprod(X,Y)),silent=TRUE)
  if (inherits(B,'try-error') || any(!is.finite(B)))
    B <- matrix(0,ncol(X),ncol(Y))
  B
}

.mdmm_init_links_structured <- function(dat, latent, links) {
  if (!dat$pD) return(links)
  X <- cbind(1,latent)
  ridge <- max(1,dat$n*.03)
  for (nm in names(dat$D)) {
    y <- dat$D[[nm]]; li <- links[[nm]]
    if (li$type %in% c('binary','count')) {
      beta <- .mdmm_opt_binary_count(X,y,rep(1,dat$n),li$beta,li$type,ridge)
      beta[-1L] <- pmax(-2,pmin(2,.5*beta[-1L]))
      li$beta <- beta
    } else if (li$type == 'nominal') {
      beta <- .mdmm_opt_nominal(X,y,rep(1,dat$n),li$beta,ridge)
      beta[,-1L] <- pmax(-2,pmin(2,.5*beta[,-1L,drop=FALSE]))
      li$beta <- beta
    } else if (li$type == 'ordinal') {
      li <- .mdmm_opt_ordinal(latent,y,rep(1,dat$n),li,ridge)
      li$beta <- pmax(-2,pmin(2,.5*li$beta))
    }
    links[[nm]] <- li
  }
  links
}

.mdmm_initialize_structured <- function(dat, pars, variant=1L) {
  if (pars$two) {
    zC <- .mdmm_init_scores(dat$C,pars$head[[1L]]$q)
    zD <- .mdmm_init_scores(.mdmm_init_discrete_embedding(dat),pars$head[[2L]]$q)
    zhead <- list(zC,zD)
    jointly <- cbind(zC/sqrt(ncol(zC)),zD/sqrt(ncol(zD)))
  } else {
    cpart <- if (dat$pC) dat$C/sqrt(dat$pC) else matrix(0,dat$n,0L)
    dpart <- .mdmm_init_discrete_embedding(dat)
    if (ncol(dpart)) dpart <- dpart/sqrt(dat$pD)
    zhead <- list(.mdmm_init_scores(cbind(cpart,dpart),pars$head[[1L]]$q))
    jointly <- zhead[[1L]]
  }
  ztail <- .mdmm_init_scores(jointly,pars$qt)
  ntail <- length(pars$tail$eta)
  tail_group <- .mdmm_init_groups(ztail,ntail,variant)
  for (a in seq_len(ntail)) {
    zz <- ztail[tail_group==a,,drop=FALSE]
    pars$tail$eta[[a]] <- colMeans(zz)
    pars$tail$psi[[a]] <- pmax(.25,colMeans(sweep(zz,2L,colMeans(zz),'-')^2)+.08)
  }
  pars$mixing[[length(pars$mixing)]] <-
    (tabulate(tail_group,nbins=ntail)+.05)/(dat$n+.05*ntail)
  head_groups <- vector("list",length(pars$head))
  for (h in seq_along(pars$head)) {
    hh <- pars$head[[h]]
    x <- zhead[[h]]
    kh <- length(hh$eta)
    grp <- .mdmm_init_groups(x,kh,variant + 3L*h)
    head_groups[[h]] <- grp
    B <- .mdmm_init_ridge(cbind(1,ztail),x)
    A <- t(B[-1L,,drop=FALSE])*.65
    for (a in seq_len(kh)) {
      idx <- which(grp==a)
      hh$A[[a]] <- A
      hh$eta[[a]] <- as.numeric(colMeans(x[idx,,drop=FALSE]) -
                                  A %*% colMeans(ztail[idx,,drop=FALSE]))
      resid <- x[idx,,drop=FALSE] -
        matrix(hh$eta[[a]],length(idx),hh$q,byrow=TRUE) -
        ztail[idx,,drop=FALSE] %*% t(A)
      hh$psi[[a]] <- pmax(.15,colMeans(resid^2)+.08)
      # Small, data-based asymmetric starts only for the skew-t family.
      # A zero empirical third moment retains a symmetric initial fit.
      if (pars$family == 'skew-t') {
        centered <- sweep(x[idx,,drop=FALSE],2L,colMeans(x[idx,,drop=FALSE]),'-')
        denom <- pmax(colMeans(centered^2),1e-5)^1.5
        hh$delta[[a]] <- as.numeric(.05*pmax(-1,pmin(1,colMeans(centered^3)/denom)))
      }
    }
    pars$mixing[[h]] <- (tabulate(grp,nbins=kh)+.05)/(dat$n+.05*kh)
    pars$head[[h]] <- hh
  }
  if (dat$pC) {
    zc <- zhead[[1L]]
    for (a in seq_along(pars$measurement$b)) {
      idx <- if (pars$two) which(head_groups[[1L]]==a) else seq_len(dat$n)
      XX <- cbind(1,zc[idx,,drop=FALSE])
      YY <- dat$C[idx,,drop=FALSE]
      B <- .mdmm_init_ridge(XX,YY)
      pars$measurement$b[[a]] <- as.numeric(B[1L,])
      pars$measurement$L[[a]] <- t(B[-1L,,drop=FALSE])
      resid <- YY - XX %*% B
      pars$measurement$psi[[a]] <- pmax(.12,colMeans(resid^2)+.05)
    }
  }
  if (dat$pD) {
    zd <- zhead[[if (pars$two) 2L else 1L]]
    pars$links <- .mdmm_init_links_structured(dat,zd,pars$links)
  }
  .mdmm_identify(pars)
}
