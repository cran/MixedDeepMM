# Input encoding; response types must be explicit to avoid silently treating
# categorical identifiers as continuous numeric data.
.mdmm_prepare <- function(data, types, schema = NULL) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (is.null(names(data)) || anyDuplicated(names(data)))
    stop("data must have unique column names.")
  if (is.null(schema)) {
    if (is.null(types) || is.null(names(types)) || anyNA(types) ||
        anyDuplicated(names(types)) ||
        !identical(sort(names(types)), sort(names(data))))
      stop("types must be a named vector specifying every column in data.")
    types <- tolower(as.character(types[names(data)]))
    if (any(!types %in% c("continuous", "binary", "ordinal", "nominal", "count")))
      stop("types may only contain continuous, binary, ordinal, nominal, count.")
    names(types) <- names(data)
    schema <- vector("list", ncol(data)); names(schema) <- names(data)
  } else {
    types <- vapply(schema, `[[`, "", "type")
    if (!setequal(names(data), names(schema)))
      stop("newdata must have the same named columns as training data.")
    data <- data[names(schema)]
  }
  C <- list(); D <- list()
  for (nm in names(types)) {
    ty <- types[[nm]]
    v <- data[[nm]]
    if (anyNA(v)) stop("Missing values in '", nm, "' are not yet supported. Use complete cases or impute explicitly.")
    if (ty == "continuous") {
      if (!is.numeric(v) || any(!is.finite(v))) stop(nm, " must be finite numeric.")
      if (is.null(schema[[nm]])) {
        center <- mean(v); scl <- stats::sd(v)
        if (!is.finite(scl) || scl < 1e-8) stop(nm, " has zero variance.")
        schema[[nm]] <- list(type = ty, center = center, scale = scl)
      }
      C[[nm]] <- (as.numeric(v) - schema[[nm]]$center) / schema[[nm]]$scale
    } else if (ty == "count") {
      if (!is.numeric(v) || any(!is.finite(v)) || any(v < 0) || any(abs(v-round(v)) > 1e-8))
        stop(nm, " must be non-negative integer-valued numeric.")
      D[[nm]] <- as.numeric(v)
      if (is.null(schema[[nm]])) schema[[nm]] <- list(type = ty)
    } else if (ty == "binary") {
      if (is.null(schema[[nm]])) {
        levels <- if (is.logical(v)) c(FALSE, TRUE) else if (is.factor(v) || is.character(v))
          unique(as.character(v)) else c(0, 1)
        if (length(levels) != 2L) stop(nm, " requires exactly two known levels.")
        schema[[nm]] <- list(type = ty, levels = as.character(levels))
      }
      y <- match(as.character(v), schema[[nm]]$levels) - 1L
      if (anyNA(y)) stop(nm, " contains a level absent from training data.")
      D[[nm]] <- y
    } else {
      if (is.null(schema[[nm]])) {
        levels <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        if (length(levels) < 2L) stop(nm, " must contain >= 2 levels.")
        if (ty == "ordinal" && !is.ordered(v))
          stop(nm, " must be an ordered factor for explicit level ordering.")
        schema[[nm]] <- list(type = ty, levels = as.character(levels))
      }
      y <- match(as.character(v), schema[[nm]]$levels)
      if (anyNA(y)) stop(nm, " contains a level absent from training data.")
      D[[nm]] <- y
    }
  }
  n <- nrow(data)
  yC <- if (length(C)) do.call(cbind, C) else matrix(numeric(0), n, 0)
  colnames(yC) <- names(C)
  list(C = yC, D = D, schema = schema, types = types, n = n, pC = ncol(yC), pD = length(D))
}

.mdmm_init_links <- function(dat, q) {
  out <- list()
  for (nm in names(dat$D)) {
    y <- dat$D[[nm]]; ty <- dat$types[[nm]]
    if (ty == "binary") {
      out[[nm]] <- list(type = ty, beta = c(qlogis(pmin(.99, pmax(.01, mean(y)))), rep(0, q)))
    } else if (ty == "count") {
      out[[nm]] <- list(type = ty, beta = c(log(pmax(.05, mean(y))), rep(0, q)))
    } else if (ty == "ordinal") {
      L <- length(dat$schema[[nm]]$levels)
      probs <- vapply(seq_len(L-1L), function(j) mean(y <= j), numeric(1))
      thresholds <- qlogis(pmin(.98, pmax(.02, probs)))
      # seq.int(2, 1) would yield c(2, 1): never iterate when there
      # is a single cutpoint (a valid two-level ordered response).
      if (length(thresholds) > 1L)
        for (j in seq.int(2L, length(thresholds)))
          thresholds[j] <- max(thresholds[j], thresholds[j-1L]+.1)
      out[[nm]] <- list(type=ty, thresholds=thresholds, beta=rep(0,q))
    } else {
      L <- length(dat$schema[[nm]]$levels)
      freqs <- tabulate(y, nbins=L) + .1
      b <- matrix(0, L-1L, q+1L)
      b[,1] <- log(freqs[seq_len(L-1L)] / freqs[L])
      out[[nm]] <- list(type=ty, beta=b)
    }
  }
  out
}

.mdmm_initialize <- function(dat, k, r, heads, family, variant=1L) {
  expected_k <- if (heads == "two") 3L else 2L
  if (!is.numeric(k) || anyNA(k) || any(!is.finite(k)) ||
      length(k) != expected_k || any(k < 1) || any(k != trunc(k)) ||
      any(k > .Machine$integer.max))
    stop("Component counts k must be positive whole numbers for the requested heads.")
  k <- as.integer(k)
  two <- heads == "two"
  if (two) {
    if (!dat$pC || !dat$pD) stop("two heads require both continuous and discrete columns.")
    if (length(k) != 3L || !is.list(r) ||
        !all(c("continuous", "discrete", "tail") %in% names(r)))
      stop("For two heads use k=c(kC,kD,kTail), r=list(continuous=qC,discrete=qD,tail=qT).")
    dims <- c(r$continuous, r$discrete, r$tail)
    if (!is.numeric(dims) || length(dims) != 3L || anyNA(dims) ||
        any(!is.finite(dims)) || any(dims != trunc(dims)))
      stop("Two-head latent dimensions must be finite whole numbers.")
    qs <- as.integer(dims[1:2]); qt <- as.integer(dims[3])
    hnames <- c("continuous", "discrete")
  } else {
    if (length(k) != 2L || is.list(r) || length(r) != 2L)
      stop("For one head use k=c(kHead,kTail), r=c(qHead,qTail).")
    if (!is.numeric(r) || anyNA(r) || any(!is.finite(r)) ||
        any(r != trunc(r)))
      stop("Single-head latent dimensions must be finite whole numbers.")
    qs <- as.integer(r[1]); qt <- as.integer(r[2])
    hnames <- "shared"
  }
  if (length(qs) != length(hnames) || anyNA(c(qs,qt)) ||
      any(c(qs,qt) <= 0) || qt >= min(qs))
    stop("Component counts and dimensions must be positive integers and tail dimension < head dimensions.")
  H <- length(qs); qx <- sum(qs)
  head <- vector("list", H)
  start <- 1L
  for (h in seq_len(H)) {
    inds <- seq.int(start, length.out=qs[h]); start <- start + qs[h]
    head[[h]] <- list(name=hnames[h], idx=inds, q=qs[h],
                      eta=lapply(seq_len(k[h]), function(a) rep(0,qs[h])),
                      A=lapply(seq_len(k[h]), function(a) matrix(0,qs[h],qt)),
                      psi=lapply(seq_len(k[h]), function(a) rep(1,qs[h])),
                      delta=lapply(seq_len(k[h]), function(a) rep(0,qs[h])))
  }
  xC <- if (two) head[[1]]$idx else head[[1]]$idx
  xD <- if (two) head[[2]]$idx else head[[1]]$idx
  LC <- matrix(0,dat$pC,length(xC))
  # M2DGMM: the observation-level continuous MFA parameters are specific
  # to the continuous-head component. M1DGMM has a shared Gaussian GLLVM
  # measurement link with conditionally independent continuous responses.
  meas <- list(b=lapply(seq_len(if (two) k[1L] else 1L),
                        function(a) rep(0,dat$pC)),
               L=lapply(seq_len(if (two) k[1L] else 1L),
                        function(a) LC),
               psi=lapply(seq_len(if (two) k[1L] else 1L),
                          function(a) rep(.7,dat$pC)))
  pars <- list(head=head, tail=list(eta=lapply(seq_len(tail(k,1)),function(a) rep(0,qt)),
                                     psi=lapply(seq_len(tail(k,1)),function(a) rep(1,qt)),
                                     delta=lapply(seq_len(tail(k,1)),function(a) rep(0,qt))),
               mixing=lapply(k,function(z)rep(1/z,z)),
               measurement=meas,
               links=.mdmm_init_links(dat,length(xD)),
               family=family, nu=8, qx=qx, qt=qt, xC=xC, xD=xD, k=k, two=two)
  .mdmm_initialize_structured(dat,pars,variant=variant)
}

.mdmm_link_prob <- function(x, link) {
  ty <- link$type
  if (ty == "binary") {
    eta <- drop(cbind(1,x) %*% link$beta)
    return(cbind(1-.mdmm_prob(eta), .mdmm_prob(eta)))
  }
  if (ty == "count") return(pmin(exp(pmin(eta <- drop(cbind(1,x) %*% link$beta), 25)),1e10))
  if (ty == "ordinal") {
    eta <- drop(x %*% link$beta)
    cp <- sapply(link$thresholds, function(a) .mdmm_prob(a - eta))
    cp <- matrix(cp, nrow=nrow(x))
    out <- cbind(cp[,1], if (ncol(cp)>1) cp[,2:ncol(cp),drop=FALSE]-cp[,1:(ncol(cp)-1L),drop=FALSE], 1-cp[,ncol(cp)])
    return(pmax(out,1e-14) / rowSums(pmax(out,1e-14)))
  }
  eta <- cbind(1,x) %*% t(link$beta)
  eta <- cbind(eta, 0)
  mx <- apply(eta,1,max)
  ex <- exp(eta-mx)
  ex / rowSums(ex)
}

.mdmm_discrete_log <- function(x, dat, i, pars) {
  if (!dat$pD) return(rep(0,nrow(x)))
  xd <- x[,pars$xD,drop=FALSE]
  z <- rep(0,nrow(x))
  for (nm in names(dat$D)) {
    link <- pars$links[[nm]]
    p <- .mdmm_link_prob(xd, link)
    yy <- dat$D[[nm]][i]
    if (link$type == "count") {
      z <- z + stats::dpois(yy, lambda=p, log=TRUE)
    } else {
      z <- z + log(pmax(p[cbind(seq_len(nrow(x)),rep(if (link$type=="binary") yy+1L else yy,nrow(x)))],1e-300))
    }
  }
  z
}
