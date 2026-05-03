test_that("point estimates are stable on reference seed", {
  set.seed(20260212)
  x <- cbind(rnorm(250), rnorm(250, sd = 1.3))
  b <- c(0.75, 0.9)
  m <- c(2, 2)
  p <- c(0.1, -0.2)

  ash <- ASH(p, x, b = b, m = m)
  lb <- LBFP(p, x, b = b)
  gl <- GLBFP(p, x, b = b, m = m)

  expect_equal(ash$estimation, 0.124444444444444, tolerance = 1e-10)
  expect_equal(lb$estimation, 0.130069621071414, tolerance = 1e-10)
  expect_equal(gl$estimation, 0.122526513628746, tolerance = 1e-10)
})

test_that("grid estimate summary is stable on reference seed", {
  set.seed(20260212)
  x <- cbind(rnorm(250), rnorm(250, sd = 1.3))
  b <- c(0.75, 0.9)
  m <- c(2, 2)

  out <- GLBFP_estimate(x, b = b, m = m, grid_size = 8)

  expect_equal(mean(out$densities), 0.0167041295031142, tolerance = 1e-10)
  expect_equal(max(out$densities), 0.105911507546916, tolerance = 1e-10)
})
