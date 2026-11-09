test_that("unbounded zero completion agrees with independent observation bases", {
  set.seed(109)
  for (d in 1:4) {
    x <- matrix(runif(20 * d, -1, 1), ncol = d)
    q <- rbind(x, rep(-20, d), rep(20, d))
    for (method in c("ASH", "LBFP", "GLBFP")) {
      b <- seq(.7, 1.2, length.out = d)
      m <- if (method == "LBFP") rep(1L, d) else rep(2L, d)
      ref <- reference_density(x, q, b, m, rep(-1, d), method)
      for (engine in c("R", "cpp_direct", "cpp_prefix", "auto")) {
        fit <- prepare_glbfp(x, b, m, method, rep(-1, d), rep(1, d), engine)
        expect_equal(predict(fit, q), ref, tolerance = 1e-11)
        expect_true(all(predict(fit, q) >= 0))
      }
    }
  }
})

test_that("ASH is literally the average of shifted coarse histograms", {
  x <- matrix(c(-.9, -.4, .1, .5, 1.4), ncol = 1)
  q <- matrix(seq(-2, 2, length.out = 31), ncol = 1)
  b <- 1.2; m <- 3L; origin <- -.9
  exact <- vapply(q[, 1], function(z) {
    mean(vapply(0:(m - 1L), function(shift) {
      a <- origin + shift * b / m
      mean(floor((x[, 1] - a) / b) == floor((z - a) / b)) / b
    }, numeric(1)))
  }, numeric(1))
  expect_equal(predict(prepare_glbfp(x, b, m, "ASH", origin, 2), q), exact, tolerance = 1e-12)
})

test_that("all compact basis densities integrate to one", {
  x <- matrix(c(0, .1, 1.2, 1.9, 2), ncol = 1)
  for (method in c("ASH", "LBFP", "GLBFP")) {
    for (mi in if (method == "LBFP") 1 else c(1, 2, 4)) {
      delta <- 1 / mi
      q <- seq(-4, 6, by = delta) + delta / 2
      y <- predict(prepare_glbfp(x, 1, mi, method, 0, 2), matrix(q))
      mass <- if (method == "ASH") sum(y) * delta else
        sum(diff(q) * (head(y, -1) + tail(y, -1)) / 2)
      expect_equal(mass, 1, tolerance = 1e-12)
      expect_equal(y[c(1, length(y))], c(0, 0))
    }
  }
})

test_that("all D_i scores match explicit deletions on the original grid", {
  x <- matrix(c(0, 0, .2, .8, 1.4, 2), ncol = 1)
  for (method in c("ASH", "LBFP", "GLBFP")) {
    mi <- if (method == "LBFP") 1 else 3
    fit <- prepare_glbfp(x, 1, mi, method, 0, 2)
    out <- compute_Di(fit)
    for (i in seq_len(nrow(x))) {
      minus <- reference_density(x[-i, , drop = FALSE], x[i, , drop = FALSE], 1, mi, 0, method)
      expect_equal(out$D[i], 1 - minus / out$density[i], tolerance = 1e-12)
      expect_equal(out$density_loo[i], minus, tolerance = 1e-12)
    }
    expect_true(all(out$D <= 1 + 1e-12))
    expect_true(all(out$D >= -1 / (nrow(x) - 1) - 1e-12))
  }
})

test_that("IID variance equals the independent sample variance of weights", {
  set.seed(119)
  x <- matrix(rnorm(100), ncol = 1)
  q <- matrix(c(-1, .1, 1))
  fit <- GLBFP_estimate(x, b = .8, m = 2, grid_points = q,
    min_vals = -5, max_vals = 5, inference = "iid")
  W <- reference_weights(x, q, .8, 2, -5)
  expect_equal(fit$sd^2, reference_iid_variance(W), tolerance = 1e-12)
  expect_equal((fit$IC[, 2] - fit$IC[, 1]) / 2, qnorm(.975) * fit$sd)
  expect_true(all(is.na(GLBFP(0, x, b = .8)$IC)))
  expect_error(GLBFP(0, x, b = .8, inference = "iid"), "explicit fixed")
  expect_error(GLBFP_estimate(x, b = .8, min_vals = -5, inference = "iid"), "fixed `max_vals`")
})

test_that("Gaussian-reference bandwidth is equivariant under changes of units", {
  set.seed(120)
  x <- cbind(rnorm(300), rnorm(300) * 2)
  scale <- c(100, .1)
  b <- compute_bi_optim(x, c(2, 3))
  expect_equal(compute_bi_optim(sweep(x, 2, scale, "*"), c(2, 3)), b * scale, tolerance = 1e-9)
  # Directly minimize equation (29), without calling G_star or equation (36).
  S <- cov(x); precision <- diag(solve(S)); d <- 2; m <- c(2, 3)
  K <- vapply(m, function(mi) integrate(function(u) sqrt((2 * mi^2 + 1 - 6*u*(1-u))/(3*mi^2)), 0, 1)$value, numeric(1))
  C0 <- sqrt(2/pi) * prod(K) * 2^(3*d/4) * pi^(d/4) * det(S)^(1/4)
  Ci <- (1/12) * (1 + 1/(2*m^2)) * 4 / sqrt(2*pi*exp(1)) * precision
  objective <- function(log_b) C0 / sqrt(nrow(x) * prod(exp(log_b))) + sum(Ci * exp(2*log_b))
  opt <- optim(log(c(1, 1)), objective, method = "BFGS", control = list(reltol = 1e-13))
  expect_equal(b, exp(opt$par), tolerance = 1e-5)
  one <- matrix(x[, 1])
  implicit <- GLBFP(0, one, m = 3)
  expect_equal(implicit$b, compute_bi_optim(one, m = 3))
})

test_that("lattice dimension is independent of observation dimension", {
  set.seed(121)
  for (N in 1:5) {
    dims <- rep(2L, N)
    a <- array(rnorm(prod(dims) * 2), c(dims, 2))
    field <- glbfp_field(a, value_dim = N + 1L, units = c("cm", "s"))
    expect_equal(field$N, N)
    expect_equal(field$d, 2L)
    expect_true(field$complete)
    b <- c(1, 1); m <- c(2, 2)
    out <- compute_Di(field, b, m, min_vals = c(-6, -6), max_vals = c(6, 6))
    ref <- compute_Di(field$values, b, m, min_vals = c(-6, -6), max_vals = c(6, 6))
    expect_equal(out$D, ref$D)
    expect_equal(out$indices, field$indices)
    expect_equal(as.matrix(as.data.frame(out)[paste0("index_", seq_len(N))]),
                 unname(field$indices), ignore_attr = TRUE)
    expect_error(GLBFP(c(0, 0), field, b, m, min_vals = c(-6, -6), max_vals = c(6, 6), inference = "iid"), "dependent-field")
  }
  a <- array(1:24, c(2, 3, 4)); a[2, 2, 2] <- NA
  expect_error(glbfp_field(a), "missing")
  field <- glbfp_field(a, na.action = "omit")
  expect_false(field$complete)
  expect_equal(nrow(field$omitted), 1L)
  expect_error(glbfp_field(1:3, rbind(c(1, 2), c(1, 2), c(2, 2))), "Duplicate")
  expect_error(glbfp_field(1:3, matrix(c(1, 2, 3.1))), "integer")
  a <- array(seq_len(2 * 3 * 4), c(2, 3, 4))
  expect_equal(glbfp_field(a, value_dim = 1)$values,
               glbfp_field(aperm(a, c(2, 3, 1)), value_dim = 3)$values)
})

test_that("grid limits do not truncate or extrapolate density support", {
  x <- matrix(c(0, .2, .9, 1.2))
  q <- matrix(seq(-2, 3, length.out = 30))
  a <- GLBFP_estimate(x, 1, 2, grid_points = q, min_vals = 0, max_vals = 1.2)
  b <- GLBFP_estimate(x, 1, 2, grid_points = q, min_vals = 0, max_vals = 10)
  expect_equal(a$densities, b$densities)
  expect_equal(GLBFP(1e300, x, b = 1)$estimation, 0)
  expect_error(GLBFP(0, x, b = 1e-300), "scale|indices")
  expect_error(GLBFP(0, x, b = 1, m = 1e20), "integer range")
  expect_error(GLBFP(1e10, matrix(rep(1e10, 3)), b = 1e-12,
    min_vals = 1e10, max_vals = 1e10 + 1), "resolve cell centers")
  expect_error(GLBFP_estimate(matrix(rnorm(40), 4, 10), b = rep(1, 10)), "million")
  expect_error(GLBFP_estimate(x, b = 1, grid_points = matrix(numeric(), 0, 1)), "at least one")
})

test_that("rectangular plots are invariant to row order", {
  set.seed(122)
  x <- cbind(rnorm(40), rnorm(40))
  g <- as.matrix(expand.grid(a = seq(-2, 2, length.out = 7), b = seq(-3, 3, length.out = 5)))
  fit <- GLBFP_estimate(x, b = c(1, 1), grid_points = g)
  shuffled <- fit; perm <- sample(nrow(g)); shuffled$grid <- g[perm, ]; shuffled$densities <- fit$densities[perm]
  p <- plotly::plotly_build(plot_density(fit, "surface"))$x$data[[1]]
  ps <- plotly::plotly_build(plot_density(shuffled, "surface"))$x$data[[1]]
  expect_equal(dim(p$z), c(5L, 7L))
  expect_equal(p$z, ps$z)
  expect_equal(p$z[2, 4], fit$densities[which(g[, 1] == p$x[4] & g[, 2] == p$y[2])])
  dup <- g; dup[1, ] <- dup[2, ]
  fit_dup <- GLBFP_estimate(x, b = c(1, 1), grid_points = dup)
  expect_false(fit_dup$is_rectangular)
  expect_error(plot_density(fit_dup, "heatmap"), "complete Cartesian")
  expect_s3_class(plot_density(fit, "heatmap"), "ggplot")
  x3 <- cbind(x, rnorm(40))
  fit3 <- GLBFP_estimate(x3, b = c(1, 1, 1), grid_size = 6)
  expect_s3_class(plot_density(fit3), "ggplot")
  expect_s3_class(plot_density(fit3, "isosurface"), "plotly")
})

test_that("prepared native objects rebuild after serialization", {
  x <- matrix(c(-1, 0, 1))
  model <- prepare_glbfp(x, 1)
  file <- tempfile(); on.exit(unlink(file))
  saveRDS(model, file)
  expect_equal(predict(readRDS(file)), predict(model))
})

test_that("tensor-product integration is normalized for d=2 and d=3", {
  set.seed(123)
  for (d in 2:3) {
    x <- matrix(runif(12 * d), ncol = d)
    for (method in c("ASH", "LBFP", "GLBFP")) {
      m <- rep(if (method == "LBFP") 1L else 2L, d)
      b <- rep(1, d); delta <- 1 / m[1]
      axis <- seq(-3, 4, by = delta) + delta/2
      grid <- as.matrix(expand.grid(rep(list(axis), d)))
      value <- predict(prepare_glbfp(x, b, m, method, rep(0,d), rep(1,d)), grid)
      expect_equal(sum(value) * delta^d, 1, tolerance = 1e-11)
    }
  }
  for (mi in c(1, 2, 10, 1e6, 1e200)) {
    numeric <- integrate(function(u) sqrt((2 + (1/mi)^2 * (1 - 6*u*(1-u)))/3), 0, 1)$value
    expect_equal(K_mi(mi), numeric, tolerance = 1e-10)
  }
})
