if (!requireNamespace("bench", quietly = TRUE)) {
  stop("Package 'bench' is required. Install it with install.packages('bench').")
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}

library(GLBFP)

set.seed(20260212)
quick <- identical(tolower(Sys.getenv("GLBFP_BENCH_QUICK")), "true")
iterations <- if (quick) 1 else 5

simulate_2d <- function(n) {
  MASS::mvrnorm(
    n = n,
    mu = c(0, 0),
    Sigma = matrix(c(1.0, 0.4, 0.4, 1.2), nrow = 2)
  )
}

run_one <- function(n, grid_size) {
  x <- simulate_2d(n)
  b <- compute_bi_optim(x, m = c(1, 1))

  bench::mark(
    ASH = ASH_estimate(x, b = b, m = c(1, 1), grid_size = grid_size),
    LBFP = LBFP_estimate(x, b = b, grid_size = grid_size),
    GLBFP = GLBFP_estimate(x, b = b, m = c(1, 1), grid_size = grid_size),
    iterations = iterations,
    check = FALSE
  )
}

settings <- if (quick) {
  expand.grid(n = 100, grid_size = 8)
} else {
  expand.grid(
    n = c(200, 500),
    grid_size = c(12, 20)
  )
}

results <- lapply(seq_len(nrow(settings)), function(i) {
  cat("Running n =", settings$n[i], "grid_size =", settings$grid_size[i], "\n")
  out <- run_one(settings$n[i], settings$grid_size[i])
  out$n <- settings$n[i]
  out$grid_size <- settings$grid_size[i]
  out
})

results <- do.call(rbind, results)
print(results)
