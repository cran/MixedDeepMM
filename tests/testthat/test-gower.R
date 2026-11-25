test_that("mixed-data Gower distance respects types and drops no rows", {
  skip_if_not_installed("cluster")
  x <- data.frame(value=c(0,1,3,4,6), yes=c(0,1,0,1,0),
                  category=factor(c("a","a","b","b","a")),
                  rating=ordered(c(1,2,1,3,2), levels=1:3))
  types <- c(value="continuous",yes="binary",category="nominal",rating="ordinal")
  d <- mixeddeepmm_gower_distance(x, types)
  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 5L)
  z <- mixeddeepmm_gower_silhouette(c(1,1,2,2,2), d)
  expect_true(is.finite(z))
  expect_equal(z, mixeddeepmm_gower_silhouette(c("A","A","B","B","B"), d))
  expect_true(is.na(mixeddeepmm_gower_silhouette(rep(1,5), d)))
  expect_error(mixeddeepmm_gower_silhouette(1:3,d), "one non-missing")
  expect_error(mixeddeepmm_gower_distance(x,types[-1]), "exactly")
})
