test_that("synthetic schema is reproducible and typed", {
  a <- simulate_mixeddeepmm(n=60,seed=1)
  b <- simulate_mixeddeepmm(n=60,seed=1)
  expect_identical(a$data,b$data)
  expect_length(a$group,60)
  expect_identical(unname(a$types),c(rep("continuous",3),rep("binary",2),
                                    "ordinal","nominal","count"))
  expect_true(is.ordered(a$data$ordinal))
  expect_error(MixedDeepMM:::.mdmm_prepare(a$data,types=a$types[-1]),"types")
  x <- a$data; x$count[1] <- -1
  expect_error(MixedDeepMM:::.mdmm_prepare(x,a$types),"non-negative")
  x <- a$data; x$ordinal <- as.character(x$ordinal)
  expect_error(MixedDeepMM:::.mdmm_prepare(x,a$types),"ordered")
})

test_that("adjusted Rand index agrees on exact matching", {
  expect_equal(mixeddeepmm_ari(c(1,1,2,2),c(9,9,3,3)),1)
  expect_true(is.finite(mixeddeepmm_ari(1:10,rep(1,10))))
})
