test_that("bundled real climate data preserve the N=2 and N=3 workflows", {
  data <- read.csv(system.file("extdata", "cangrd_quebec_1991_2020.csv", package = "GLBFP"))
  expect_equal(nrow(data), 2760L)
  expect_equal(sort(unique(data$year)), 1991:2020)
  expect_equal(nrow(unique(data[c("i", "j")])), 92L)
  expect_false(anyDuplicated(data[c("i", "j", "year")]) > 0L)
  expect_true(all(is.finite(as.matrix(data[c("tmin_anomaly_c", "tmax_anomaly_c")]))))
  fields <- lapply(c(FALSE, TRUE), function(with_time) {
    selected <- if (with_time) data else data[data$year == 2020L, ]
    glbfp_field(as.matrix(selected[c("tmin_anomaly_c", "tmax_anomaly_c")]),
      indices = as.matrix(selected[if (with_time) c("i", "j", "year") else c("i", "j")]),
      units = c("deg C", "deg C"), dependence = "dependent")
  })
  expect_equal(vapply(fields, function(x) x$N, integer(1)), c(2L, 3L))
  expect_equal(vapply(fields, function(x) x$d, integer(1)), c(2L, 2L))
  expect_equal(vapply(fields, function(x) nrow(x$values), integer(1)), c(92L, 2760L))
  expect_true(all(!vapply(fields, function(x) x$complete, logical(1))))
  for (field in fields) {
    fit <- GLBFP_estimate(field, b = c(.6, .6), m = c(3, 3), grid_size = 5,
      min_vals = c(-8, -8), max_vals = c(8, 8))
    ref <- GLBFP_estimate(field$values, b = c(.6, .6), m = c(3, 3),
      grid_points = fit$grid, min_vals = c(-8, -8), max_vals = c(8, 8), engine = "R")
    expect_equal(fit$densities, ref$densities, tolerance = 1e-12)
    expect_true(all(is.na(fit$sd)))
    expect_error(GLBFP(c(0, 0), field, b = c(.6, .6), min_vals = c(-8, -8),
      max_vals = c(8, 8), inference = "iid"), "dependent-field")
  }
})
