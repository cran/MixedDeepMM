test_that("MR is invariant to relabelling and detects errors", {
  truth <- factor(c("A", "A", "B", "B", "C", "C"))
  expect_equal(mixeddeepmm_mr(truth, c(9, 9, 4, 4, 7, 7)), 0)
  expect_equal(mixeddeepmm_mr(truth, c(4, 4, 7, 7, 9, 9)), 0)
  expect_equal(mixeddeepmm_mr(c("A", "A", "B", "B"), c(9, 9, 4, 9)), 0.25)
  expect_equal(mixeddeepmm_mr("A", "other"), 0)
})

test_that("MR handles unequal class and cluster counts", {
  expect_equal(mixeddeepmm_mr(c("A", "A", "B", "B"), c(1, 2, 3, 3)), 0.25)
  expect_equal(mixeddeepmm_mr(c("A", "A", "B", "B"), c(1, 1, 1, 1)), 0.5)
})

test_that("MR optimises assignment globally rather than class by class", {
  # Rows: A=(4,3), B=(3,0). Optimal assignment: A->2, B->1 = 6/10.
  truth <- c(rep("A", 7), rep("B", 3))
  predicted <- c(rep(1, 4), rep(2, 3), rep(1, 3))
  expect_equal(mixeddeepmm_mr(truth, predicted), 0.4)
})

test_that("MR rejects missing or incompatible inputs", {
  expect_error(mixeddeepmm_mr(1:2, 1), "equal-length")
  expect_error(mixeddeepmm_mr(c(1, NA), c(1, 2)), "non-missing")
  expect_error(mixeddeepmm_mr(integer(), integer()), "nonempty")
})
