test_that("discrete_uniform_moments returns correct mean and variance", {
  a <- 1; b <- 6
  m <- discrete_uniform_moments(a, b)
  expect_equal(m$mean, mean(1:6))
  expect_equal(m$variance, var(1:6) * (5/6)) # population correction
})