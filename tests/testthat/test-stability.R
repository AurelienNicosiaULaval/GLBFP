test_that("custom irregular grids are supported", {
  set.seed(123)
  data <- cbind(rnorm(120), rnorm(120))
  b <- c(0.8, 0.8)
  m <- c(1, 1)

  grid_points <- cbind(
    seq(-2, 2, length.out = 30),
    seq(-2, 2, length.out = 30)^2 / 4
  )

  fit <- ASH_estimate(data, b = b, m = m, grid_points = grid_points)
  expect_false(fit$is_rectangular)
  expect_equal(nrow(fit$grid), 30)

  p_contour <- plot(fit, contour = TRUE)
  p_surface <- plot(fit, contour = FALSE)

  expect_s3_class(p_contour, "ggplot")
  expect_true(inherits(p_surface, "plotly"))
})

test_that("grid metadata is consistent for regular grids", {
  set.seed(9)
  data <- cbind(rnorm(80), rnorm(80))
  fit <- LBFP_estimate(data, b = c(0.9, 0.9), grid_size = 10)

  expect_true(fit$is_rectangular)
  expect_equal(fit$grid_dims, c(10L, 10L))
  expect_equal(nrow(fit$grid), 100)
})
