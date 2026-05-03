test_that("functions reject missing values", {
  x <- cbind(rnorm(30), rnorm(30))
  x[1, 1] <- NA_real_

  expect_error(compute_bi_optim(x, m = c(1, 1)), "missing values")
  expect_error(LBFP(c(0, 0), x, b = c(1, 1)), "missing values")
  expect_error(GLBFP(c(0, 0), x, b = c(1, 1), m = c(1, 1)), "missing values")
  expect_error(ASH(c(0, 0), x, b = c(1, 1), m = c(1, 1)), "missing values")
})

test_that("functions reject non-numeric data with clear errors", {
  x <- data.frame(a = rnorm(30), b = letters[seq_len(30)])

  expect_error(LBFP(c(0, 0), x, b = c(1, 1)), "numeric columns")
  expect_error(GLBFP_estimate(x, b = c(1, 1), m = c(1, 1), grid_size = 5), "numeric columns")
})

test_that("compute_bi_optim handles degenerate covariance with ridge stabilization", {
  x <- cbind(rep(1, 20), rep(1, 20))
  b <- compute_bi_optim(x, m = c(1, 1))

  expect_equal(length(b), 2)
  expect_true(all(is.finite(b)))
  expect_true(all(b > 0))
})

test_that("constant data require explicit non-degenerate estimation bounds", {
  x <- matrix(rep(1, 20), ncol = 1)

  expect_error(
    GLBFP(1, x, b = 0.5, m = 1),
    "constant data require explicit non-degenerate bounds"
  )

  fit <- GLBFP(1, x, b = 0.5, m = 1, min_vals = 0, max_vals = 2)
  expect_s3_class(fit, "GLBFP")
  expect_true(is.finite(fit$estimation))
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
