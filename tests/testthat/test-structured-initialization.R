test_that("mixed-response initialization is data-driven and never consumes RNG", {
  sim <- simulate_mixeddeepmm(n=60L,seed=51L)
  dat <- MixedDeepMM:::.mdmm_prepare(sim$data,sim$types)
  set.seed(1234)
  before <- .Random.seed
  p1 <- MixedDeepMM:::.mdmm_initialize(dat,c(2L,2L),c(2L,1L),
                                        "single","gaussian",variant=3L)
  expect_identical(.Random.seed,before)
  set.seed(77)
  other <- .Random.seed
  p2 <- MixedDeepMM:::.mdmm_initialize(dat,c(2L,2L),c(2L,1L),
                                        "single","gaussian",variant=3L)
  expect_identical(.Random.seed,other)
  expect_equal(p1,p2,tolerance=1e-12)
  expect_true(any(abs(unlist(p1$measurement$L))>1e-8))
  expect_true(any(abs(unlist(lapply(p1$head,function(z)z$A)))>1e-8))
  expect_true(all(unlist(lapply(p1$head,function(z)z$psi))>0))
  expect_true(all(unlist(p1$mixing)>0))
  expect_true(any(abs(unlist(lapply(p1$links,function(z)z$beta)))>1e-8))
})

test_that("two-head deterministic initialization has correct shapes", {
  sim <- simulate_mixeddeepmm(n=50L,seed=52L)
  dat <- MixedDeepMM:::.mdmm_prepare(sim$data,sim$types)
  p <- MixedDeepMM:::.mdmm_initialize(dat,c(2L,2L,2L),
      list(continuous=2L,discrete=2L,tail=1L),"two","t",variant=5L)
  expect_equal(length(p$head),2L)
  expect_length(p$measurement$L,2L)
  expect_equal(dim(p$measurement$L[[1L]]),c(dat$pC,2L))
  expect_equal(dim(p$measurement$L[[2L]]),c(dat$pC,2L))
  expect_length(p$tail$psi,2L)
  expect_equal(length(p$links$b1$beta),3L)
  expect_equal(vapply(p$mixing,sum,numeric(1)),rep(1,3),tolerance=1e-10)
  expect_equal(p$nu,8)
  expect_true(all(is.finite(unlist(p$measurement))))
})

test_that("mixed starts report a deterministic initialization variant", {
  sim <- simulate_mixeddeepmm(n=35L,seed=53L)
  fit <- mixeddeepmm(sim$data,sim$types,"gaussian","single",
    layers=list(k=c(2L,2L),r=c(2L,1L)),
    control=list(max_iter=2L,min_iter=2L,mc_samples=6L,
                 mc_eval=10L,nstart=2L,seed=10L))
  expect_identical(fit$init_method,"structured_data_driven")
  expect_equal(fit$starts$variant,c(11L,12L))
  expect_true(is.finite(fit$loglik))
})
