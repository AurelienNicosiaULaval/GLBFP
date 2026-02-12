`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || !nzchar(x)) y else x

bench_dir <- dirname(sys.frame(1)$ofile %||% "")
if (!nzchar(bench_dir) || !file.exists(file.path(bench_dir, "sim_scenarios.R"))) {
  pkg_bench_dir <- system.file("bench", package = "GLBFP")
  if (nzchar(pkg_bench_dir) && file.exists(file.path(pkg_bench_dir, "sim_scenarios.R"))) {
    bench_dir <- pkg_bench_dir
  }
}
if (!nzchar(bench_dir) || !file.exists(file.path(bench_dir, "sim_scenarios.R"))) {
  bench_dir <- file.path("inst", "bench")
}
source(file.path(bench_dir, "sim_scenarios.R"))

run_validation_suite <- function(reps = 3, seed = 20260212, grid_size_1d = 80, grid_size_2d = 20) {
  set.seed(seed)

  methods <- c("ASH", "LBFP", "GLBFP", "KDE")
  out <- list()
  row_id <- 1L

  run_one_1d <- function(scenario, n) {
    x <- simulate_1d(scenario, n)
    x_mat <- matrix(x, ncol = 1)
    b <- compute_bi_optim(x_mat, m = 1)

    t_ash <- system.time(ash <- ASH_estimate(x_mat, b = b, m = 1, grid_size = grid_size_1d))["elapsed"]
    t_lbfp <- system.time(lbfp <- LBFP_estimate(x_mat, b = b, grid_size = grid_size_1d))["elapsed"]
    t_glbfp <- system.time(glbfp <- GLBFP_estimate(x_mat, b = b, m = 1, grid_size = grid_size_1d))["elapsed"]
    t_kde <- system.time(kde <- stats::density(x, n = grid_size_1d))["elapsed"]

    data.frame(
      scenario = scenario,
      dimension = 1,
      method = methods,
      ise = c(
        ise_1d(ash$grid[, 1], ash$densities, scenario),
        ise_1d(lbfp$grid[, 1], lbfp$densities, scenario),
        ise_1d(glbfp$grid[, 1], glbfp$densities, scenario),
        ise_1d(kde$x, kde$y, scenario)
      ),
      elapsed = as.numeric(c(t_ash, t_lbfp, t_glbfp, t_kde))
    )
  }

  run_one_2d <- function(scenario, n) {
    x <- simulate_2d(scenario, n)
    b <- compute_bi_optim(x, m = c(1, 1))

    t_ash <- system.time(ash <- ASH_estimate(x, b = b, m = c(1, 1), grid_size = grid_size_2d))["elapsed"]
    t_lbfp <- system.time(lbfp <- LBFP_estimate(x, b = b, grid_size = grid_size_2d))["elapsed"]
    t_glbfp <- system.time(glbfp <- GLBFP_estimate(x, b = b, m = c(1, 1), grid_size = grid_size_2d))["elapsed"]

    lims <- c(range(x[, 1]), range(x[, 2]))
    t_kde <- system.time(kde <- MASS::kde2d(x[, 1], x[, 2], n = grid_size_2d, lims = lims))["elapsed"]

    kde_grid <- expand.grid(kde$x, kde$y)

    data.frame(
      scenario = scenario,
      dimension = 2,
      method = methods,
      ise = c(
        ise_2d(ash$grid[, 1], ash$grid[, 2], ash$densities, scenario),
        ise_2d(lbfp$grid[, 1], lbfp$grid[, 2], lbfp$densities, scenario),
        ise_2d(glbfp$grid[, 1], glbfp$grid[, 2], glbfp$densities, scenario),
        ise_2d(kde_grid[, 1], kde_grid[, 2], as.vector(kde$z), scenario)
      ),
      elapsed = as.numeric(c(t_ash, t_lbfp, t_glbfp, t_kde))
    )
  }

  for (rep_idx in seq_len(reps)) {
    for (scenario in c("S1", "S2", "S3")) {
      out[[row_id]] <- transform(run_one_1d(scenario, n = 200), rep = rep_idx)
      row_id <- row_id + 1L
    }
    for (scenario in c("S4", "S5")) {
      out[[row_id]] <- transform(run_one_2d(scenario, n = 250), rep = rep_idx)
      row_id <- row_id + 1L
    }
  }

  metrics <- do.call(rbind, out)

  sensitivity <- local({
    x <- simulate_2d("S4", 250)
    b <- compute_bi_optim(x, m = c(1, 1))

    run_cfg <- function(b_scale, m_val) {
      fit <- GLBFP_estimate(x, b = b * b_scale, m = c(m_val, m_val), grid_size = grid_size_2d)
      c(
        b_scale = b_scale,
        m = m_val,
        mean_density = mean(fit$densities),
        max_density = max(fit$densities)
      )
    }

    as.data.frame(rbind(
      run_cfg(0.8, 1),
      run_cfg(1.0, 1),
      run_cfg(1.2, 1),
      run_cfg(1.0, 2)
    ))
  })

  list(metrics = metrics, sensitivity = sensitivity)
}

summarise_validation <- function(result) {
  mean_tbl <- aggregate(
    cbind(ise, elapsed) ~ scenario + dimension + method,
    data = result$metrics,
    FUN = mean
  )
  median_tbl <- aggregate(
    cbind(ise, elapsed) ~ scenario + dimension + method,
    data = result$metrics,
    FUN = stats::median
  )

  names(mean_tbl)[names(mean_tbl) == "ise"] <- "ise_mean"
  names(mean_tbl)[names(mean_tbl) == "elapsed"] <- "elapsed_mean"
  names(median_tbl)[names(median_tbl) == "ise"] <- "ise_median"
  names(median_tbl)[names(median_tbl) == "elapsed"] <- "elapsed_median"

  merge(
    mean_tbl,
    median_tbl,
    by = c("scenario", "dimension", "method"),
    sort = TRUE
  )
}

if (sys.nframe() == 0) {
  result <- run_validation_suite()
  summary_tbl <- summarise_validation(result)

  print(summary_tbl)
  cat("\nSensitivity summary\n")
  print(result$sensitivity)
}
