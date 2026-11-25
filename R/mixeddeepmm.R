#' Fit a continuous or mixed-response deep mixture model.
#'
#' `types` assigns one of continuous, binary, ordinal, nominal or count to each
#' named response. For continuous data only, uses the supplied DGMM/RDMM/DStMM
#' algorithms. For mixed data, fits a new conditional-measurement latent-mixture
#' model with Monte Carlo EM. In the two-head Gaussian case continuous
#' observation MFA loadings/intercepts/variances depend on the continuous-head
#' component. The shared-tail layer is a Gaussian mixture with freely estimated
#' diagonal covariance, which marginalizes the last factor layer for a shallow
#' tail. This is paper-aligned, not the full published MDGMM estimation algorithm.
#' Student-t and skew-t mixed kernels extend the Gaussian specification.
#'
#' @param data Data frame or numeric matrix.
#' @param types Named character vector with a type for each column, or NULL for
#'   entirely numeric continuous data.
#' @param family "gaussian", "t", or "skew-t".
#' @param heads "single" or "two". Two heads are for mixed responses only.
#' @param layers List containing k and r. For mixed single-head use
#'   `list(k=c(2,2),r=c(2,1))`; for two heads use
#'   `list(k=c(2,2,2),r=list(continuous=2,discrete=2,tail=1))`.
#' @param control Named list of algorithm settings. See README.
#' @return An object of class mixeddeepmm.
#' @importFrom stats cov cutree dist factanal filter hclust kmeans optimize plogis
#' @importFrom stats predict princomp qlogis rbinom rgamma rnorm rpois uniroot var
#' @importFrom utils head tail
#' @export
mixeddeepmm <- function(data, types=NULL, family=c("gaussian","t","skew-t"),
                        heads=c("single","two"), layers, control=list()) {
  family <- match.arg(family)
  heads <- match.arg(heads)
  if (missing(layers) || !is.list(layers) || !all(c("k","r") %in% names(layers)))
    stop("layers must be a list with k and r.")
  defaults <- list(max_iter=35L,min_iter=8L,tol=0.002,mc_samples=48L,
                   mc_eval=192L,nstart=1L,seed=123L,verbose=FALSE,
                   psi_floor=.03,ridge=.001,nu=8,estimate_nu=TRUE,
                   init="kmeans",init_est="factanal")
  if (!is.list(control) || (length(control) && is.null(names(control))))
    stop("control must be a named list.")
  # Convenient aliases for the manuscript's EM naming conventions.
  if (!is.null(control$it)) { control$max_iter <- control$it; control$it <- NULL }
  if (!is.null(control$eps)) { control$tol <- control$eps; control$eps <- NULL }
  if (!is.null(control$M)) { control$mc_samples <- control$M; control$M <- NULL }
  unknown <- setdiff(names(control),names(defaults))
  if (length(unknown)) stop("Unknown control option(s): ",paste(unknown,collapse=", "))
  co <- utils::modifyList(defaults,control)
  numeric_controls <- c("max_iter","min_iter","tol","mc_samples",
                        "mc_eval","nstart","psi_floor","ridge","nu")
  if (!all(vapply(co[numeric_controls],function(v)
      is.numeric(v) && length(v)==1L && is.finite(v),logical(1))))
    stop("Numerical control settings must be finite numeric scalars.")
  integer_controls <- c("max_iter","min_iter","mc_samples","mc_eval","nstart")
  if (!all(vapply(co[integer_controls],function(v)
      v==trunc(v) && v<=.Machine$integer.max,logical(1))))
    stop("Iteration counts, MC sizes and nstart must be positive integers.")
  if (co$max_iter < 1 || co$min_iter < 1 || co$mc_samples < 2 ||
      co$mc_eval < 2 || co$nstart < 1 || co$psi_floor <= 0 ||
      co$ridge < 0 || co$nu <= 2 || co$tol <= 0)
    stop("Invalid control values (nu must be >2, MC sizes >=2, floors/tolerance >0).")
  if (!is.null(co$seed) && (!is.numeric(co$seed) || length(co$seed)!=1L ||
      !is.finite(co$seed) || co$seed != trunc(co$seed) ||
      abs(co$seed) + co$nstart - 1 > .Machine$integer.max))
    stop("seed must be an integer within the R integer range, or NULL.")
  if (!all(vapply(co[c("estimate_nu","verbose")],function(v)
      is.logical(v) && length(v)==1L && !is.na(v),logical(1))))
    stop("estimate_nu and verbose must be TRUE or FALSE.")
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (is.null(types)) {
    if (!all(vapply(data,is.numeric,logical(1))))
      stop("Supply named types for any non-numeric/categorical responses.")
    types <- stats::setNames(rep("continuous",ncol(data)),names(data))
  }
  if (is.null(names(types)) || anyNA(types) ||
      anyDuplicated(names(types)) || !setequal(names(types),names(data)))
    stop("types must explicitly name every column in data without missing values.")
  if (all(tolower(types)=="continuous")) {
    if (heads!="single") stop("Continuous-only models use heads='single'.")
    return(.mdmm_fit_continuous(data,types,family,layers,co,match.call()))
  }
  dat <- .mdmm_prepare(data,types)
  if (family=="skew-t" && !requireNamespace("GIGrvg",quietly=TRUE))
    stop("Mixed skew-t fitting requires optional package 'GIGrvg'.")
  fits <- vector("list",as.integer(co$nstart)); starts <- data.frame()
  for (s in seq_len(as.integer(co$nstart))) {
    if (!is.null(co$seed)) set.seed(as.integer(co$seed)+s-1L)
    # Distinct deterministic structured starts (no random parameter draws).
    # The seed independently governs MCEM sampling; external scripts changing
    # the seed also explore alternative data-driven initial partitions.
    variant <- if (is.null(co$seed)) s else
      as.integer((abs(as.double(co$seed)+s-1L) %% 99991) + 1L)
    ans <- try(.mdmm_fit_one(dat,layers,heads,family,co,variant=variant),silent=TRUE)
    if (inherits(ans,"try-error")) {
      starts <- rbind(starts,data.frame(start=s,variant=variant,logLik=NA_real_,error=as.character(ans)[1]))
    } else {
      fits[[s]] <- ans
      starts <- rbind(starts,data.frame(start=s,variant=variant,logLik=ans$loglik,error=""))
    }
  }
  valid <- which(vapply(fits,function(x)!is.null(x) && is.finite(x$loglik),logical(1)))
  if (!length(valid)) stop("All mixed-response starts failed: ",paste(starts$error,collapse="; "))
  best <- fits[[valid[which.max(vapply(fits[valid],`[[`,numeric(1),"loglik"))]]]
  best$starts <- starts
  best$call <- match.call()
  best$control <- co
  class(best) <- "mixeddeepmm"
  best
}

.mdmm_fit_one <- function(dat,layers,heads,family,co,variant=1L) {
  pars <- .mdmm_initialize(dat,layers$k,layers$r,heads,family,variant=variant)
  pars$nu <- co$nu
  logtrace <- numeric(0)
  for (iter in seq_len(as.integer(co$max_iter))) {
    est <- .mdmm_estep(dat,pars,as.integer(co$mc_samples))
    logtrace <- c(logtrace,est$loglik)
    if (!is.finite(est$loglik)) stop("Nonfinite Monte Carlo log likelihood.")
    if (isTRUE(co$verbose)) message(sprintf("mixed %s iter %d: logLik=%.3f",
                                            family,iter,est$loglik))
    nextpars <- .mdmm_update(dat,pars,est,co$ridge,co$psi_floor,co$estimate_nu)
    # Parameter updates are MCEM, so ascent of sampled likelihood is NOT guaranteed.
    pars <- nextpars
    if (iter >= max(2L,co$min_iter) && length(logtrace)>=6L) {
      old <- mean(tail(logtrace,6)[1:3]); new <- mean(tail(logtrace,3))
      if (abs(new-old)/(abs(old)+1) < co$tol) break
    }
  }
  # Independent Monte Carlo likelihood estimate for BIC/model comparison.
  check <- .mdmm_estep(dat,pars,as.integer(co$mc_eval),save_draws=FALSE)
  free <- .mdmm_parameters(dat,pars,co$estimate_nu)
  LL <- check$loglik
  entropy <- -sum(check$tau*log(pmax(check$tau,1e-300)))
  # This is the paper-motivated evaluation layer: M1 common first layer;
  # M2 shared-tail layer. Do not use external labels to choose this layer.
  paper_layer <- if (heads=="two") ncol(check$paths) else 1L
  paper_prob <- matrix(0,dat$n,pars$k[paper_layer])
  for (a in seq_len(ncol(paper_prob)))
    paper_prob[,a] <- rowSums(check$tau[,
                        check$paths[,paper_layer]==a,drop=FALSE])
  list(engine="mixed_mcem", family=family,heads=heads,layers=layers,
       init_method="structured_data_driven", init_variant=variant,
       parameterization="paper_aligned_shallow_mixed",
       params=pars, schema=dat$schema,
       data=dat,posterior=check$tau,paths=check$paths,first_layer=check$first_layer,
       paper_layer=paper_layer,paper_posterior=paper_prob,
       cluster=max.col(paper_prob,ties.method="first"),
       path_cluster=max.col(check$tau,ties.method="first"),
       loglik=LL,mc_se=check$mc_se,parameters=free,
       bic=-2*LL+free*log(dat$n),
       aic=-2*LL+2*free,
       icl=-2*LL+free*log(dat$n)+2*entropy,
       entropy=entropy,trace=logtrace,
       iterations=length(logtrace),
       warning=if (family!="gaussian")
         "Monte Carlo mixed-kernel approximation; examine numerical integration and model sensitivity."
         else "Conditional Gaussian/GLLVM head with mixture tail; MCEM approximation, not the original NSEP/estimation algorithm.")
}

.mdmm_fit_continuous <- function(data,types,family,layers,co,call) {
  dat <- .mdmm_prepare(data,types)
  y <- dat$C
  k <- as.integer(layers$k); r <- as.integer(layers$r)
  if (length(k)<1L || length(k)>3L || length(k)!=length(r) ||
      any(r <= 0) || any(diff(c(ncol(y),r)) >= 0))
    stop("Continuous layers require 1-3 dimensions strictly decreasing from observed p.")
  if (family=="gaussian" && (!requireNamespace("mvtnorm",quietly=TRUE) ||
      !requireNamespace("corpcor",quietly=TRUE)))
    stop("Continuous Gaussian legacy engine needs packages 'mvtnorm' and 'corpcor'.")
  if (family=="skew-t" && !requireNamespace("GIGrvg",quietly=TRUE))
    stop("Continuous skew-t engine requires optional 'GIGrvg'.")
  starts <- data.frame(); fits <- vector("list",as.integer(co$nstart))
  for (s in seq_len(as.integer(co$nstart))) {
    seed <- if (is.null(co$seed)) NULL else as.integer(co$seed)+s-1L
    z <- try(switch(family,
      gaussian=deepgmm(y,layers=length(k),k=k,r=r,it=co$max_iter,
                      eps=co$tol,init=co$init,init_est=co$init_est,
                      seed=seed,scale=FALSE,psi_floor=co$psi_floor),
      t=robustdeepgmm(y,layers=length(k),k=k,r=r,it=co$max_iter,
                       eps=co$tol,init=co$init,init_est=co$init_est,
                       seed=seed,scale=FALSE,nu=co$nu,estimate_nu=co$estimate_nu,
                       method="ecm",psi_floor=co$psi_floor,
                       min_iter=co$min_iter),
      `skew-t`=dstmm(y,layers=length(k),k=k,r=r,it=co$max_iter,
                     eps=co$tol,init=co$init,init_est=co$init_est,
                     seed=seed,scale=FALSE,nu=co$nu,estimate_nu=co$estimate_nu,
                     M=co$mc_samples,psi_floor=co$psi_floor,
                     min_iter=co$min_iter)),silent=TRUE)
    if (inherits(z,"try-error")) starts <- rbind(starts,data.frame(start=s,logLik=NA_real_,error=as.character(z)[1]))
    else {
      fits[[s]] <- z
      ll <- if (!is.null(z$loglik)) as.numeric(z$loglik) else tail(z$lik,1)
      starts <- rbind(starts,data.frame(start=s,logLik=ll,error=""))
    }
  }
  valid <- which(is.finite(starts$logLik))
  if (!length(valid)) stop("All continuous starts failed: ",paste(starts$error,collapse="; "))
  best <- valid[which.max(starts$logLik[valid])]
  z <- fits[[best]]; ll <- starts$logLik[best]
  npar <- z$h
  b <- if (length(npar)==1L && is.finite(npar)) -2*ll+npar*log(nrow(y)) else z$bic
  out <- list(engine="legacy_continuous", model=z, family=family,
              heads="single",layers=layers,
              schema=dat$schema, data=dat, loglik=ll,parameters=npar,
              bic=b,aic=-2*ll+2*npar,icl=if (!is.null(z$icl_bic)) z$icl_bic else NA_real_,
              mc_se=NA_real_,trace=z$lik,
              iterations=if (!is.null(z$iterations)) z$iterations else length(z$lik),
              starts=starts,control=co,call=call)
  class(out) <- "mixeddeepmm"
  out
}
