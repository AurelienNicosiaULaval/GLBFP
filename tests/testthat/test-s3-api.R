test_that("point fit objects keep backward-compatible classes", {
  set.seed(2024)
  data <- cbind(rnorm(150), rnorm(150))
  fit <- GLBFP(c(0, 0), data, b = c(0.8, 0.8), m = c(1, 1))

  expect_true(inherits(fit, "glbfp_fit"))
  expect_true(inherits(fit, "GLBFP"))

  s <- summary(fit)
  expect_s3_class(s, "summary.glbfp_fit")
  expect_equal(predict(fit), fit$estimation)
  expect_error(predict(fit, newdata = data[1:2, ]), "not supported")
})

test_that("grid fit objects expose summary and nearest-grid prediction", {
  set.seed(2025)
  data <- cbind(rnorm(120), rnorm(120))
  fit <- GLBFP_estimate(data, b = c(0.9, 0.9), m = c(1, 1), grid_size = 9)

  expect_true(inherits(fit, "glbfp_grid"))
  expect_true(inherits(fit, "GLBFP_estimate"))

  s <- summary(fit)
  expect_s3_class(s, "summary.glbfp_grid")

  pred_all <- predict(fit)
  expect_equal(length(pred_all), nrow(fit$grid))

  newdata <- fit$grid[1:3, , drop = FALSE]
  pred_new <- predict(fit, newdata = newdata)
  expect_equal(length(pred_new), 3)
  expect_true(all(is.finite(pred_new)))
})
