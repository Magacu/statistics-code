test_that("discrete_uniform PMF and CDF work", {
  a <- 1; b <- 6
  expect_equal(sum(discrete_uniform_pmf(1:6, a, b)), 1)
  expect_equal(discrete_uniform_cdf(3, a, b), 3/6)
})