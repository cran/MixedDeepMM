test_that("single head mixed Gaussian MCEM creates normalized posteriors", {
  sim <- simulate_mixeddeepmm(n=40,seed=10)
  fit <- mixeddeepmm(sim$data,sim$types,"gaussian","single",
    layers=list(k=c(2L,2L),r=c(2L,1L)),
    control=list(max_iter=2,min_iter=2,mc_samples=8,mc_eval=12,seed=1))
  expect_s3_class(fit,"mixeddeepmm")
  pp <- predict(fit,type="posterior")
  expect_equal(dim(pp),c(40L,4L))
  expect_equal(rowSums(pp),rep(1,40),tolerance=1e-6)
  expect_equal(dim(predict(fit,type="layer",layer=1L)),c(40L,2L))
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_false("ess" %in% names(fit))
  expect_false(any(grepl("ess",names(mixeddeepmm_diagnostics(fit)),fixed=TRUE)))
})

test_that("two-head structure shares a tail", {
  sim <- simulate_mixeddeepmm(n=40,seed=22)
  fit <- mixeddeepmm(sim$data,sim$types,"gaussian","two",
    layers=list(k=c(2L,2L,2L),r=list(continuous=2L,discrete=2L,tail=1L)),
    control=list(max_iter=2,min_iter=2,mc_samples=8,mc_eval=12,seed=2))
  expect_equal(dim(fit$posterior),c(40L,8L))
  expect_true(is.finite(fit$bic))
  expect_equal(rowSums(fit$posterior),rep(1,40),tolerance=1e-6)
})

test_that("continuous-only measurement collapses to exact observed density", {
  set.seed(25)
  dat <- data.frame(a=rnorm(30),b=rnorm(30),c=rnorm(30))
  tps <- c(a="continuous",b="continuous",c="continuous")
  d <- MixedDeepMM:::.mdmm_prepare(dat,tps)
  p <- MixedDeepMM:::.mdmm_initialize(d,c(2L,2L),c(2L,1L),"single","gaussian")
  ee <- MixedDeepMM:::.mdmm_estep(d,p,M=5,save_draws=FALSE)
  paths <- MixedDeepMM:::.mdmm_paths(p$k)
  direct <- 0
  for(i in seq_len(30)) {
    lp <- vapply(seq_len(nrow(paths)),function(s) {
      cc <- MixedDeepMM:::.mdmm_path_cache(p,paths[s,])
      MixedDeepMM:::.mdmm_cont_posterior(d$C[i,],cc,"gaussian",8)$logd+cc$logpi
    },numeric(1))
    direct <- direct+MixedDeepMM:::.mdmm_logsumexp(lp)
  }
  expect_equal(ee$loglik,direct,tolerance=1e-8)
  expect_equal(ee$mc_se,0,tolerance=1e-12)
})

test_that("mixed t fits return finite log likelihood", {
  sim <- simulate_mixeddeepmm(n=40,seed=11)
  fit <- mixeddeepmm(sim$data,sim$types,"t","single",
    layers=list(k=c(2L,2L),r=c(2L,1L)),
    control=list(max_iter=2,min_iter=2,mc_samples=8,mc_eval=12,seed=1))
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$family, "t")
})

test_that("mixed skew-t fits when GIGrvg is available", {
  skip_if_not_installed("GIGrvg")
  sim <- simulate_mixeddeepmm(n=40,seed=9)
  fit <- mixeddeepmm(sim$data,sim$types,"skew-t","single",
    layers=list(k=c(2L,2L),r=c(2L,1L)),
    control=list(max_iter=2,min_iter=2,mc_samples=8,mc_eval=12,seed=1))
  expect_true(is.finite(fit$loglik))
})
