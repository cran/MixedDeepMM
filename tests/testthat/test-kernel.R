test_that("Gaussian and Student-t kernels follow known analytic formulas", {
  dat <- simulate_mixeddeepmm(n=40,seed=12)
  d <- MixedDeepMM:::.mdmm_prepare(dat$data,dat$types)
  p <- MixedDeepMM:::.mdmm_initialize(d,c(2L,2L),c(2L,1L),"single","gaussian")
  cc <- MixedDeepMM:::.mdmm_path_cache(p,c(1L,1L))
  y <- d$C[1,]
  got <- MixedDeepMM:::.mdmm_cont_posterior(y,cc,"gaussian",nu=8)
  ref <- MixedDeepMM:::.mdmm_mvnorm_log(y,cc$meanC,cc$sigma)
  expect_equal(got$logd,ref,tolerance=1e-8)
  expect_true(is.finite(MixedDeepMM:::.mdmm_cont_posterior(y,cc,"t",8)$logd))
  expect_equal(MixedDeepMM:::.mdmm_cont_posterior(y,cc,"skew-t",8)$logd,
               MixedDeepMM:::.mdmm_cont_posterior(y,cc,"t",8)$logd,
               tolerance=1e-8) # zero skew delta
})

test_that("anchor transformations preserve Gaussian path and link predictions", {
  sim <- simulate_mixeddeepmm(n=50,seed=8)
  d <- MixedDeepMM:::.mdmm_prepare(sim$data,sim$types)
  p <- MixedDeepMM:::.mdmm_initialize(d,c(2L,2L),c(2L,1L),"single","gaussian")
  b <- lapply(seq_len(4),function(s) MixedDeepMM:::.mdmm_path_cache(p,MixedDeepMM:::.mdmm_paths(p$k)[s,]))
  p2 <- MixedDeepMM:::.mdmm_identify(p)
  a <- lapply(seq_len(4),function(s) MixedDeepMM:::.mdmm_path_cache(p2,MixedDeepMM:::.mdmm_paths(p2$k)[s,]))
  for(s in seq_len(4)) {
    expect_equal(b[[s]]$meanC,a[[s]]$meanC,tolerance=1e-9)
    expect_equal(b[[s]]$sigma,a[[s]]$sigma,tolerance=1e-9)
  }
})

test_that("Gaussian and t likelihoods agree as nu grows", {
  set.seed(21)
  sim <- simulate_mixeddeepmm(n=42,seed=21)
  d <- MixedDeepMM:::.mdmm_prepare(sim$data,sim$types)
  p <- MixedDeepMM:::.mdmm_initialize(d,c(2L,2L),c(2L,1L),"single","gaussian")
  cc <- MixedDeepMM:::.mdmm_path_cache(p,c(2L,1L))
  y <- d$C[2,]
  a <- MixedDeepMM:::.mdmm_cont_posterior(y,cc,"gaussian",8)$logd
  b <- MixedDeepMM:::.mdmm_cont_posterior(y,cc,"t",1e7)$logd
  expect_equal(a,b,tolerance=1e-5)
})
