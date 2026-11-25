#' Simulate a mixed response dataset with known clusters.
#'
#' Generates three continuous, two binary, one ordinal, one nominal and one
#' count response from three groups sharing a two-dimensional latent vector.
#' Class labels are provided separately, not as a fitted response.
#'
#' @param n Number of independent observations.
#' @param seed Random seed.
#' @return List with `data`, `types`, and generating `group`.
#' @export
simulate_mixeddeepmm <- function(n=600L,seed=2026L) {
  if (length(n)!=1 || n<30 || n!=as.integer(n)) stop("n must be integer >=30.")
  set.seed(seed)
  group <- sample.int(3,n,replace=TRUE,prob=c(.35,.4,.25))
  cent <- rbind(c(-1.1,.5),c(1,.9),c(.1,-1.2))
  latent <- cent[group,,drop=FALSE]+matrix(rnorm(n*2,sd=.75),n,2)
  X <- latent[,1]; Z <- latent[,2]
  logistic <- function(x) plogis(x)
  c1 <- 1.0*X+.35*Z+rnorm(n,sd=.55)
  c2 <- -.4*X+1.0*Z+rnorm(n,sd=.65)
  c3 <- .7*X-.8*Z+rnorm(n,sd=.6)
  b1 <- rbinom(n,1,logistic(-.1+1.1*X-.45*Z))
  b2 <- rbinom(n,1,logistic(.3-.5*X+.9*Z))
  ordinal_lp <- .65*X+.75*Z+stats::rlogis(n)
  ord <- cut(ordinal_lp,breaks=c(-Inf,-.6,.65,Inf),
             labels=c("low","middle","high"),ordered_result=TRUE)
  eta <- cbind(.4+.8*X,.3-.9*Z,rep(0,n))
  ee <- exp(eta-apply(eta,1,max)); pr <- ee/rowSums(ee)
  nominal <- factor(apply(pr,1,function(p)sample.int(3,1,prob=p)),
                    levels=1:3,labels=c("A","B","C"))
  count <- rpois(n,lambda=exp(.2+.3*X-.4*Z))
  data <- data.frame(c1,c2,c3,b1,b2,ordinal=ord,nominal=nominal,count=count)
  types <- c(c1="continuous",c2="continuous",c3="continuous",
             b1="binary",b2="binary",ordinal="ordinal",
             nominal="nominal",count="count")
  list(data=data,types=types,group=group,latent=latent)
}
