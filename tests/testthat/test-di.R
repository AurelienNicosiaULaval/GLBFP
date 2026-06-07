test_that("compute_Di matches explicit GLBFP leave-one-out recomputation", {
  set.seed(20260607)
  x <- cbind(rnorm(50), rnorm(50, sd = 1.2))
  b <- c(0.8, 0.9)
  m <- c(2, 2)
  min_vals <- apply(x, 2, min)
  max_vals <- apply(x, 2, max)

  out <- compute_Di(x, b = b, m = m, estimator = "GLBFP",
                    min_vals = min_vals, max_vals = max_vals)

  expect_s3_class(out, "glbfp_di")
  expect_equal(length(out$D), nrow(x))
  expect_true(all(is.finite(out$D)))

  for (i in seq_len(6L)) {
    f_minus <- GLBFP(
      x = x[i, ],
      data = x[-i, , drop = FALSE],
      b = b,
      m = m,
      min_vals = min_vals,
      max_vals = max_vals
    )$estimation
    expect_equal(out$D[i], 1 - f_minus / out$density[i], tolerance = 1e-12)
  }
})

test_that("compute_Di supports LBFP and ASH", {
  set.seed(20260608)
  x <- cbind(rnorm(40), rnorm(40))
  b <- c(0.7, 0.7)
  min_vals <- apply(x, 2, min)
  max_vals <- apply(x, 2, max)

  lbfp <- compute_Di(x, b = b, estimator = "LBFP",
                     min_vals = min_vals, max_vals = max_vals)
  ash <- compute_Di(x, b = b, m = c(2, 2), estimator = "ASH",
                    min_vals = min_vals, max_vals = max_vals)

  expect_s3_class(lbfp, "glbfp_di")
  expect_s3_class(ash, "glbfp_di")
  expect_true(all(is.finite(lbfp$D)))
  expect_true(all(is.finite(ash$D)))

  s <- summary(lbfp)
  expect_s3_class(s, "summary.glbfp_di")
})
