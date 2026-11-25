test_that("continuous degrees-of-freedom is free only when estimated", {
  k <- c(3L, 2L)
  r_full <- c(13L, 4L, 2L)
  fixed_rd <- MixedDeepMM:::.rdmm_parameter_count(k, r_full, "common", FALSE)
  fitted_rd <- MixedDeepMM:::.rdmm_parameter_count(k, r_full, "common", TRUE)
  expect_equal(fitted_rd - fixed_rd, 1L)
  fixed_st <- MixedDeepMM:::.dstmm_parameter_count(k, r_full, 1L, "common", FALSE)
  fitted_st <- MixedDeepMM:::.dstmm_parameter_count(k, r_full, 1L, "common", TRUE)
  expect_equal(fitted_st - fixed_st, 1L)
  expect_equal(MixedDeepMM:::.rdmm_parameter_count(k, r_full, "first_layer", TRUE) - fixed_rd, 3L)
  expect_equal(MixedDeepMM:::.rdmm_parameter_count(k, r_full, "pathway", TRUE) - fixed_rd, 6L)
})

test_that("conditional scale moments can move nu away from its initial value", {
  tau <- matrix(1/2, nrow=20L, ncol=2L)
  paths <- matrix(c(1L, 1L, 2L, 1L), ncol=2L, byrow=TRUE)
  rd_estep <- list(paths=paths, tau=tau,
                   elogw=matrix(-0.5, nrow=20L, ncol=2L),
                   ew=matrix(1, nrow=20L, ncol=2L))
  rd_nu <- MixedDeepMM:::.rdmm_update_nu(rd_estep, rep(8,2), "common", c(2L,1L), c(2.05,200))
  expect_true(all(is.finite(rd_nu)))
  expect_true(all(rd_nu > 2))
  expect_false(isTRUE(all.equal(rd_nu, rep(8,2))))
  st_estep <- list(paths=paths, tau=tau,
                   elogH=matrix(0.5, nrow=20L, ncol=2L),
                   einvh=matrix(1, nrow=20L, ncol=2L))
  st_nu <- MixedDeepMM:::.dstmm_update_nu(st_estep, rep(8,2), "common", c(2L,1L), c(2.05,200))
  expect_true(all(is.finite(st_nu)))
  expect_true(all(st_nu > 2))
  expect_false(isTRUE(all.equal(st_nu, rep(8,2))))
})
