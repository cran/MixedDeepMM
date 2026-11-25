test_that("two-head continuous measurement is component specific", {
  sim <- simulate_mixeddeepmm(n=45L,seed=72L)
  dat <- MixedDeepMM:::.mdmm_prepare(sim$data,sim$types)
  pars <- MixedDeepMM:::.mdmm_initialize(dat,c(2L,2L,2L),
    list(continuous=2L,discrete=2L,tail=1L),"two","gaussian",variant=9L)
  expect_length(pars$measurement$b,2L)
  expect_length(pars$measurement$L,2L)
  expect_length(pars$measurement$psi,2L)
  expect_equal(dim(pars$measurement$L[[1]]),c(dat$pC,2L))
  expect_equal(dim(pars$measurement$L[[2]]),c(dat$pC,2L))
  expect_true(all(unlist(pars$tail$psi)>0))
  paths <- MixedDeepMM:::.mdmm_paths(pars$k)
  s <- paths[1L,]
  cache1 <- MixedDeepMM:::.mdmm_path_cache(pars,s)
  pars2 <- pars
  pars2$measurement$b[[1]][1] <- pars2$measurement$b[[1]][1]+3
  cache2 <- MixedDeepMM:::.mdmm_path_cache(pars2,s)
  expect_equal(unname(cache2$meanC[1]-cache1$meanC[1]),3,tolerance=1e-8)
  # Change to continuous component 2: parameter of component 1 must not matter.
  s[1] <- 2L
  cache3 <- MixedDeepMM:::.mdmm_path_cache(pars,s)
  cache4 <- MixedDeepMM:::.mdmm_path_cache(pars2,s)
  expect_equal(cache3$meanC,cache4$meanC,tolerance=1e-10)
})

test_that("paper clustering level changes with architecture, explicit layers work", {
  sim <- simulate_mixeddeepmm(n=35L,seed=73L)
  ctl <- list(max_iter=2L,min_iter=2L,mc_samples=6L,mc_eval=9L,seed=3L)
  single <- mixeddeepmm(sim$data,sim$types,"gaussian","single",
    layers=list(k=c(2L,2L),r=c(2L,1L)),control=ctl)
  two <- mixeddeepmm(sim$data,sim$types,"gaussian","two",
    layers=list(k=c(2L,2L,2L),
      r=list(continuous=2L,discrete=2L,tail=1L)),control=ctl)
  expect_equal(predict(single,type="cluster"),predict(single,type="cluster",layer=1L))
  expect_equal(predict(two,type="cluster"),predict(two,type="cluster",layer=3L))
  expect_equal(two$cluster, predict(two,type="cluster"))
  expect_equal(predict(two,type="cluster",cluster_level="first"),
               predict(two,type="cluster",layer=1L))
  expect_equal(predict(two,type="layer",cluster_level="paper"),
               predict(two,type="layer",layer=3L))
  expect_equal(predict(two,type="layer",cluster_level="shared_tail"),
               predict(two,type="layer",layer=3L))
  expect_error(predict(single,type="cluster",cluster_level="discrete_head"),"two-head")
  expect_false("ess" %in% names(two))
  expect_false("ess_trace" %in% names(two))
  expect_false(any(grepl("ess",names(mixeddeepmm_diagnostics(two)),fixed=TRUE)))
})

test_that("tail covariance changes head and measurement covariance", {
  sim <- simulate_mixeddeepmm(n=36L,seed=74L)
  dat <- MixedDeepMM:::.mdmm_prepare(sim$data,sim$types)
  p <- MixedDeepMM:::.mdmm_initialize(dat,c(2L,2L,2L),
    list(continuous=2L,discrete=2L,tail=1L),"two","gaussian")
  s <- MixedDeepMM:::.mdmm_paths(p$k)[1L,]
  a <- MixedDeepMM:::.mdmm_path_cache(p,s)
  p$tail$psi[[1]] <- p$tail$psi[[1]]*2
  b <- MixedDeepMM:::.mdmm_path_cache(p,s)
  expect_true(any(abs(a$sigma-b$sigma)>1e-8))
})
