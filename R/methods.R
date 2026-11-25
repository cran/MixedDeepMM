#' @export
print.mixeddeepmm <- function(x, ...) {
  cat("MixedDeepMM", if(x$engine=="mixed_mcem")"(MCEM)" else "(continuous legacy)","\n")
  cat("Family:",x$family," | Heads:",x$heads,"\n")
  if (!is.null(x$init_method))
    cat("Initialization:",x$init_method," | Variant:",x$init_variant,"\n")
  cat("LogLik:",format(round(x$loglik,3))," | BIC:",format(round(x$bic,3)),"\n")
  if (is.finite(x$mc_se)) cat("Independent integration MC SE (logLik):",signif(x$mc_se,3),"\n")
  invisible(x)
}

#' @export
summary.mixeddeepmm <- function(object, ...) {
  x <- list(family=object$family,heads=object$heads,layers=object$layers,
            engine=object$engine,loglik=object$loglik,bic=object$bic,
            icl=object$icl,npar=object$parameters,mc_se=object$mc_se,
            iterations=object$iterations,
            starts=object$starts,
            init_method=object$init_method,init_variant=object$init_variant)
  class(x) <- "summary.mixeddeepmm"
  x
}

#' @export
print.summary.mixeddeepmm <- function(x,...) {
  cat("MixedDeepMM fit summary\n")
  cat(sprintf("Family: %s; heads: %s; engine: %s\n",x$family,x$heads,x$engine))
  cat(sprintf("logLik: %.4f; BIC: %.4f; parameters (nominal): %d\n",x$loglik,x$bic,x$npar))
  if (is.finite(x$mc_se)) cat(sprintf("Monte Carlo SE(logLik): %.4f\n",x$mc_se))
  cat("Iterations:",x$iterations,"\n")
  if (!is.null(x$init_method)) cat("Initialization:",x$init_method,
                                   "variant",x$init_variant,"\n")
  print(x$starts,row.names=FALSE)
  invisible(x)
}

#' @export
logLik.mixeddeepmm <- function(object, ...) {
  structure(object$loglik, class="logLik",df=object$parameters,
            nobs=object$data$n)
}

.mdmm_get_posterior <- function(object,dat,mc_samples=NULL) {
  if (object$engine=="mixed_mcem") {
    M <- .mdmm_prediction_mc(object,mc_samples)
    return(.mdmm_estep(dat,object$params,M,save_draws=FALSE)$tau)
  }
  z <- object$model
  y <- dat$C
  if (object$family=="gaussian") {
    ee <- compute.lik(y,nrow(y),z$k,z$mu,z$H,z$psi,z$w)
    return(ee$ps.y)
  }
  if (object$family=="t") return(predict.rdmm(z,newdata=y,type="posterior"))
  predict.dstmm(z,newdata=y,type="posterior")
}

.mdmm_link_fitted <- function(object,dat,M) {
  est <- .mdmm_estep(dat,object$params,M,save_draws=TRUE)
  out <- list()
  for (nm in names(dat$D)) {
    link <- object$params$links[[nm]]
    levels <- if (link$type %in% c("binary","ordinal","nominal")) {
      if (link$type=="binary") dat$schema[[nm]]$levels else dat$schema[[nm]]$levels
    } else NULL
    ncol <- if (link$type=="count") 1L else length(levels)
    pp <- matrix(0,dat$n,ncol)
    for (rec in est$records) {
      p <- .mdmm_link_prob(rec$z[,object$params$xD,drop=FALSE],link)
      if (link$type=="count") pp[rec$i,1] <- pp[rec$i,1]+sum(rec$w*p)
      else pp[rec$i,] <- pp[rec$i,]+colSums(p*rec$w)
    }
    if (!is.null(levels)) colnames(pp) <- levels
    out[[nm]] <- if(link$type=="count") as.vector(pp[,1]) else pp
  }
  out
}

# Validate the draw count before invoking conditional Monte Carlo prediction.
.mdmm_prediction_mc <- function(object, mc_samples=NULL) {
  M <- if (is.null(mc_samples)) object$control$mc_eval else mc_samples
  if (!is.numeric(M) || length(M) != 1L || is.na(M) ||
      !is.finite(M) || M < 2 || M != trunc(M) || M > .Machine$integer.max)
    stop("mc_samples must be one integer >= 2.")
  as.integer(M)
}

#' @export
predict.mixeddeepmm <- function(object, newdata=NULL,
                                type=c("cluster","posterior","path","layer","response"),
                                layer=NULL,mc_samples=NULL,
                                cluster_level=c("paper","first","shared_tail","continuous_head","discrete_head"), ...) {
  type <- match.arg(type)
  dat <- if (is.null(newdata)) object$data else .mdmm_prepare(newdata,types=NULL,schema=object$schema)
  if (type=="response") {
    if (object$engine!="mixed_mcem") stop("Response probabilities apply to mixed response models.")
    M <- .mdmm_prediction_mc(object,mc_samples)
    return(.mdmm_link_fitted(object,dat,M))
  }
  if (object$engine=="mixed_mcem" && !is.null(mc_samples))
    .mdmm_prediction_mc(object,mc_samples)
  prob <- if(is.null(newdata) && object$engine=="mixed_mcem" && is.null(mc_samples))
    object$posterior else .mdmm_get_posterior(object,dat,mc_samples)
  if (type=="posterior") return(prob)
  if (object$engine=="mixed_mcem") paths <- object$paths else {
    paths <- .rdmm_paths(as.integer(object$layers$k))
  }
  if (type=="path") return(max.col(prob,ties.method="first"))
  if (is.null(layer)) {
    cluster_level <- match.arg(cluster_level)
    layer <- switch(cluster_level,
      paper=if(object$heads=="two") ncol(paths) else 1L,
      first=1L,
      shared_tail=ncol(paths),
      continuous_head=1L,
      discrete_head=if(object$heads=="two") 2L else
        stop("discrete_head requires a two-head model."))
  }
  if (!is.numeric(layer) || length(layer)!=1L || is.na(layer) ||
      layer!=trunc(layer) || layer<1 || layer>ncol(paths)) stop("Invalid layer index.")
  nclass <- max(paths[,layer])
  marginal <- matrix(0,nrow(prob),nclass)
  for (a in seq_len(nclass)) marginal[,a] <- rowSums(prob[,paths[,layer]==a,drop=FALSE])
  if (type=="layer") return(marginal)
  max.col(marginal,ties.method="first")
}

#' @export
fitted.mixeddeepmm <- function(object,...) {
  if (object$engine=="mixed_mcem") return(predict(object,type="response",...))
  predict(object,type="cluster",...)
}

#' @export
plot.mixeddeepmm <- function(x, ...) {
  ll <- as.numeric(x$trace)
  if (!length(ll)) stop("No likelihood trace to plot.")
  graphics::plot(seq_along(ll),ll,type="b",xlab="Iteration",
                 ylab="Estimated log likelihood",main=paste("MixedDeepMM",x$family),...)
  invisible(x)
}

#' Inspect likelihood history and Monte Carlo precision.
#' @export
mixeddeepmm_diagnostics <- function(object) {
  stopifnot(inherits(object,"mixeddeepmm"))
  list(family=object$family,heads=object$heads,
       iterations=object$iterations,
       loglik=object$loglik,
       monte_carlo_se_loglik=object$mc_se,
       likelihood_trace=object$trace,starts=object$starts,
       init_method=object$init_method,init_variant=object$init_variant,
       warning=if(!is.null(object$warning))object$warning else
         "Continuous-source algorithm; review posterior and numerical settings.")
}

#' Adjusted Rand index for external cluster labels.
#' @export
mixeddeepmm_ari <- function(truth, predicted) {
  if (length(truth)!=length(predicted) || anyNA(truth) || anyNA(predicted))
    stop("truth and predicted need equal non-missing lengths.")
  tab <- table(truth,predicted)
  choose2 <- function(x) x*(x-1)/2
  a <- sum(choose2(tab)); b <- sum(choose2(rowSums(tab)))
  c <- sum(choose2(colSums(tab))); total <- choose2(sum(tab))
  if (total==0) return(1)
  expect <- b*c/total; denom <- (b+c)/2-expect
  if (abs(denom)<1e-12) return(if(a==b && b==c)1 else 0)
  (a-expect)/denom
}
