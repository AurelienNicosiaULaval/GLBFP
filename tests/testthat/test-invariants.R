trapz <- function(x, y) {
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

test_that("grid estimators return non-negative densities", {
  set.seed(42)
  x <- cbind(rnorm(300), rnorm(300))
  b <- c(0.7, 0.7)
  m <- c(1, 1)

  out_ash <- ASH_estimate(x, b = b, m = m, grid_size = 12)
  out_lbfp <- LBFP_estimate(x, b = b, grid_size = 12)
  out_glbfp <- GLBFP_estimate(x, b = b, m = m, grid_size = 12)

  expect_true(all(out_ash$densities >= 0))
  expect_true(all(out_lbfp$densities >= 0))
  expect_true(all(out_glbfp$densities >= 0))

  expect_equal(length(out_ash$densities), nrow(out_ash$grid))
  expect_equal(length(out_lbfp$densities), nrow(out_lbfp$grid))
  expect_equal(length(out_glbfp$densities), nrow(out_glbfp$grid))
})

test_that("1D grid estimates integrate close to one", {
  set.seed(2026)
  x <- matrix(rnorm(350), ncol = 1)
  b <- compute_bi_optim(x, m = 1)

  out_lbfp <- LBFP_estimate(x, b = b, grid_size = 120)
  out_glbfp <- GLBFP_estimate(x, b = b, m = 1, grid_size = 120)

  int_lbfp <- trapz(out_lbfp$grid[, 1], out_lbfp$densities)
  int_glbfp <- trapz(out_glbfp$grid[, 1], out_glbfp$densities)

  expect_true(abs(int_lbfp - 1) < 0.25)
  expect_true(abs(int_glbfp - 1) < 0.35)
})
