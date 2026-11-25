test_that("mixed fits omit the end-of-run convergence status", {
  sim <- simulate_mixeddeepmm(n=40L, seed=31L)
  fit <- mixeddeepmm(sim$data,sim$types,family="gaussian",heads="single",
    layers=list(k=c(2L,2L),r=c(2L,1L)),
    control=list(max_iter=2L,min_iter=2L,mc_samples=6L,mc_eval=10L,seed=1L))
  expect_false("converged" %in% names(fit))
  expect_false("converged" %in% names(summary(fit)))
  expect_false("converged" %in% names(mixeddeepmm_diagnostics(fit)))
  expect_false(any(grepl("converg",capture.output(print(summary(fit))),ignore.case=TRUE)))
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$iterations,length(fit$trace))
})

test_that("two-level ordinal response works in single and two-head mixed fits", {
  sim <- simulate_mixeddeepmm(n=40L,seed=33L)
  sim$data$ordinal <- ordered(ifelse(sim$data$ordinal=="low","low","other"),
                              levels=c("low","other"))
  d <- MixedDeepMM:::.mdmm_prepare(sim$data,sim$types)
  pars <- MixedDeepMM:::.mdmm_initialize(d,c(2L,2L),c(2L,1L),"single","gaussian")
  expect_length(pars$links$ordinal$thresholds,1L)
  expect_true(all(is.finite(pars$links$ordinal$thresholds)))
  co <- list(max_iter=2L,min_iter=2L,mc_samples=6L,mc_eval=10L,seed=9L)
  single <- mixeddeepmm(sim$data,sim$types,"gaussian","single",
                        layers=list(k=c(2L,2L),r=c(2L,1L)),control=co)
  two <- mixeddeepmm(sim$data,sim$types,"gaussian","two",
                     layers=list(k=c(2L,2L,2L),
                                 r=list(continuous=2L,discrete=2L,tail=1L)),
                     control=co)
  for (fit in list(single,two)) {
    expect_true(is.finite(fit$loglik))
    pp <- predict(fit,type="posterior")
    expect_equal(rowSums(pp),rep(1,40),tolerance=1e-6)
    probs <- predict(fit,type="response",mc_samples=8L)
    for (nm in c("b1","b2","ordinal","nominal")) {
      expect_true(all(is.finite(probs[[nm]])))
      expect_equal(rowSums(probs[[nm]]),rep(1,40),tolerance=1e-6)
    }
    expect_equal(ncol(probs$ordinal),2L)
    expect_true(all(is.finite(probs$count) & probs$count >= 0))
    permuted <- sim$data[rev(names(sim$data))]
    newpp <- predict(fit,newdata=permuted,type="posterior",mc_samples=8L)
    expect_equal(dim(newpp),dim(pp))
    expect_equal(rowSums(newpp),rep(1,40),tolerance=1e-6)
    expect_error(predict(fit,type="posterior",mc_samples=1L),"mc_samples")
    expect_error(predict(fit,type="response",mc_samples=NA_real_),"mc_samples")
  }
})

test_that("mixed single-head model also supports entirely discrete responses", {
  sim <- simulate_mixeddeepmm(n=40L,seed=34L)
  cols <- c("b1","b2","ordinal","nominal","count")
  fit <- mixeddeepmm(sim$data[cols],sim$types[cols],"gaussian","single",
                     layers=list(k=c(2L,2L),r=c(2L,1L)),
                     control=list(max_iter=2L,min_iter=2L,
                                  mc_samples=8L,mc_eval=12L,seed=4L))
  expect_true(is.finite(fit$loglik))
  expect_equal(rowSums(fit$posterior),rep(1,40),tolerance=1e-6)
  expect_true(all(is.finite(predict(fit,type="response",mc_samples=8L)$count)))
})

test_that("mixed dimensions are checked before integer coercion", {
  sim <- simulate_mixeddeepmm(n=40L,seed=35L)
  co <- list(max_iter=2L,min_iter=2L,mc_samples=6L,mc_eval=8L)
  expect_error(mixeddeepmm(sim$data,sim$types,layers=list(k=c(2.5,2),r=c(2,1)),
                            control=co),"whole numbers")
  expect_error(mixeddeepmm(sim$data,sim$types,layers=list(k=c(2,2),r=c(2.5,1)),
                            control=co),"whole numbers")
})
