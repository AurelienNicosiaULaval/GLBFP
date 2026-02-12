# Reproducible benchmark scenarios for GLBFP.

simulate_1d <- function(scenario, n) {
  switch(
    scenario,
    S1 = rnorm(n, mean = 0, sd = 1),
    S2 = ifelse(runif(n) < 0.5, rnorm(n, -2, 1), rnorm(n, 2, 1)),
    S3 = rt(n, df = 4),
    stop("Unknown 1D scenario")
  )
}

true_density_1d <- function(x, scenario) {
  switch(
    scenario,
    S1 = dnorm(x, mean = 0, sd = 1),
    S2 = 0.5 * dnorm(x, -2, 1) + 0.5 * dnorm(x, 2, 1),
    S3 = dt(x, df = 4),
    stop("Unknown 1D scenario")
  )
}

dbivnorm <- function(x, y, mu, Sigma) {
  Sigma_inv <- solve(Sigma)
  Sigma_det <- det(Sigma)
  diff1 <- x - mu[1]
  diff2 <- y - mu[2]
  quad <- Sigma_inv[1, 1] * diff1^2 +
    2 * Sigma_inv[1, 2] * diff1 * diff2 +
    Sigma_inv[2, 2] * diff2^2
  (1 / (2 * pi * sqrt(Sigma_det))) * exp(-0.5 * quad)
}

simulate_2d <- function(scenario, n) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for 2D simulation scenarios.")
  }

  if (scenario == "S4") {
    z <- runif(n) < 0.5
    x <- matrix(NA_real_, nrow = n, ncol = 2)
    x[z, ] <- MASS::mvrnorm(sum(z), mu = c(-1, -1), Sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2))
    x[!z, ] <- MASS::mvrnorm(sum(!z), mu = c(1.5, 1.5), Sigma = matrix(c(1, -0.4, -0.4, 1), 2, 2))
    return(x)
  }

  if (scenario == "S5") {
    return(MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1.2, 0.8, 0.8, 1), 2, 2)))
  }

  stop("Unknown 2D scenario")
}

true_density_2d <- function(x, y, scenario) {
  if (scenario == "S4") {
    d1 <- dbivnorm(x, y, mu = c(-1, -1), Sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2))
    d2 <- dbivnorm(x, y, mu = c(1.5, 1.5), Sigma = matrix(c(1, -0.4, -0.4, 1), 2, 2))
    return(0.5 * d1 + 0.5 * d2)
  }

  if (scenario == "S5") {
    return(dbivnorm(x, y, mu = c(0, 0), Sigma = matrix(c(1.2, 0.8, 0.8, 1), 2, 2)))
  }

  stop("Unknown 2D scenario")
}

trapz <- function(x, y) {
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

ise_1d <- function(grid, est, scenario) {
  f0 <- true_density_1d(grid, scenario)
  trapz(grid, (est - f0)^2)
}

ise_2d <- function(grid_x, grid_y, est, scenario) {
  f0 <- true_density_2d(grid_x, grid_y, scenario)
  dx <- median(diff(sort(unique(grid_x))))
  dy <- median(diff(sort(unique(grid_y))))
  sum((est - f0)^2) * dx * dy
}
