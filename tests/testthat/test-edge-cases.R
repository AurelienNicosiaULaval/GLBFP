test_that("functions reject missing values", {
  x <- cbind(rnorm(30), rnorm(30))
  x[1, 1] <- NA_real_

  expect_error(compute_bi_optim(x, m = c(1, 1)), "missing values")
  expect_error(LBFP(c(0, 0), x, b = c(1, 1)), "missing values")
  expect_error(GLBFP(c(0, 0), x, b = c(1, 1), m = c(1, 1)), "missing values")
  expect_error(ASH(c(0, 0), x, b = c(1, 1), m = c(1, 1)), "missing values")
})

test_that("compute_bi_optim handles degenerate covariance with ridge stabilization", {
  x <- cbind(rep(1, 20), rep(1, 20))
  b <- compute_bi_optim(x, m = c(1, 1))

  expect_equal(length(b), 2)
  expect_true(all(is.finite(b)))
  expect_true(all(b > 0))
})

test_that("estimators work on tiny samples and outliers", {
  small <- matrix(c(-1, 0, 1, 0, 1, -1), ncol = 2)
  outlier <- rbind(small, c(1000, -1000))

  b <- c(1, 1)

  expect_s3_class(LBFP(c(0, 0), small, b = b), "LBFP")
  expect_s3_class(GLBFP(c(0, 0), small, b = b, m = c(1, 1)), "GLBFP")
  expect_s3_class(ASH(c(0, 0), small, b = b, m = c(1, 1)), "ASH")

  fit <- GLBFP(c(0, 0), outlier, b = c(10, 10), m = c(1, 1))
  expect_true(is.finite(fit$estimation))
})
