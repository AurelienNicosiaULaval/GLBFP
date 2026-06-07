#' Summarize GLBFP fit objects
#'
#' Summarizes objects returned by [ASH()], [LBFP()], [GLBFP()] and
#' their grid counterparts.
#'
#' @param object A fitted object of class `"glbfp_fit"` or `"glbfp_grid"`.
#' @param ... Additional arguments (unused).
#'
#' @return A list with class `"summary.glbfp_fit"` or `"summary.glbfp_grid"`.
#' @export
summary.glbfp_fit <- function(object, ...) {
  out <- list(
    method = object$method,
    dimension = object$dimension,
    point = object$x,
    estimation = object$estimation,
    sd = object$sd %||% NA_real_,
    IC = object$IC %||% c(NA_real_, NA_real_),
    b = object$b,
    m = object$m %||% NULL,
    u = object$u %||% NULL,
    cell_index = object$cell_index %||% NULL,
    visited = object$visited %||% NA_integer_,
    prefix_nodes = object$prefix_nodes %||% NA_integer_
  )
  class(out) <- "summary.glbfp_fit"
  out
}

#' @export
print.summary.glbfp_fit <- function(x, ...) {
  cat("Method:", x$method, "\n")
  cat("Dimension:", x$dimension, "\n")
  cat("Point:", paste(x$point, collapse = ", "), "\n")
  cat("Estimation:", x$estimation, "\n")
  if (!is.na(x$sd)) {
    cat("Standard error:", x$sd, "\n")
  }
  if (all(is.finite(x$IC))) {
    cat("95% CI:", paste(x$IC, collapse = ", "), "\n")
  }
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  if (!is.null(x$m)) {
    cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
  }
  if (!is.null(x$u)) {
    cat("Relative grid coordinate (u):", paste(x$u, collapse = ", "), "\n")
  }
  if (is.finite(x$visited)) {
    cat("Visited cells:", x$visited, "\n")
  }
  if (is.finite(x$prefix_nodes)) {
    cat("Prefix nodes:", x$prefix_nodes, "\n")
  }
  invisible(x)
}

#' @export
summary.glbfp_grid <- function(object, ...) {
  out <- list(
    method = object$method,
    dimension = ncol(object$grid),
    n_grid = nrow(object$grid),
    grid_dims = object$grid_dims,
    is_rectangular = object$is_rectangular,
    b = object$b,
    m = object$m %||% NULL,
    density_min = min(object$densities),
    density_q25 = unname(stats::quantile(object$densities, 0.25)),
    density_median = stats::median(object$densities),
    density_mean = mean(object$densities),
    density_q75 = unname(stats::quantile(object$densities, 0.75)),
    density_max = max(object$densities),
    density_zero = sum(object$densities == 0),
    sd_min = if (!is.null(object$sd)) min(object$sd, na.rm = TRUE) else NA_real_,
    sd_median = if (!is.null(object$sd)) stats::median(object$sd, na.rm = TRUE) else NA_real_,
    sd_max = if (!is.null(object$sd)) max(object$sd, na.rm = TRUE) else NA_real_,
    visited_median = if (!is.null(object$visited)) stats::median(object$visited, na.rm = TRUE) else NA_real_,
    prefix_nodes_median = if (!is.null(object$prefix_nodes)) stats::median(object$prefix_nodes, na.rm = TRUE) else NA_real_
  )
  class(out) <- "summary.glbfp_grid"
  out
}

#' @export
print.summary.glbfp_grid <- function(x, ...) {
  cat("Method:", x$method, "\n")
  cat("Dimension:", x$dimension, "\n")
  cat("Grid points:", x$n_grid, "\n")
  cat("Grid type:", if (isTRUE(x$is_rectangular)) "rectangular" else "irregular", "\n")
  cat("Grid dimensions:", paste(x$grid_dims, collapse = " x "), "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  if (!is.null(x$m)) {
    cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
  }
  cat("Density range:", paste(c(x$density_min, x$density_max), collapse = " to "), "\n")
  cat("Density quartiles:", paste(c(x$density_q25, x$density_median, x$density_q75), collapse = ", "), "\n")
  cat("Density median:", x$density_median, "\n")
  cat("Density mean:", x$density_mean, "\n")
  cat("Zero densities:", x$density_zero, "\n")
  if (is.finite(x$sd_median)) {
    cat("Standard error median:", x$sd_median, "\n")
  }
  if (is.finite(x$visited_median)) {
    cat("Median visited cells:", x$visited_median, "\n")
  }
  if (is.finite(x$prefix_nodes_median)) {
    cat("Median prefix nodes:", x$prefix_nodes_median, "\n")
  }
  invisible(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
