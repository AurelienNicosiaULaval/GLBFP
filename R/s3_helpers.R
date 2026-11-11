# Shared print and plot helpers.

#' @keywords internal
glbfp_format_number <- function(x, digits = 6L) {
  format(signif(x, digits), trim = TRUE)
}

#' @keywords internal
glbfp_format_vector <- function(x, digits = 6L) {
  paste(glbfp_format_number(x, digits = digits), collapse = ", ")
}

#' @keywords internal
glbfp_print_fit <- function(x, label) {
  cat(label, " Density Estimation:\n", sep = "")
  cat("Point:", paste0("(", glbfp_format_vector(x$x), ")"), "\n")
  cat("Estimated density:", glbfp_format_number(x$estimation), "\n")
  if (!is.null(x$sd) && is.finite(x$sd)) {
    cat("Estimated standard error:", glbfp_format_number(x$sd), "\n")
  }
  if (!is.null(x$IC) && all(is.finite(x$IC))) {
    label <- if (identical(x$inference, "stationary_mixing")) "95% asymptotic interval for density:" else "95% normal interval for expected estimate:"
    cat(label, glbfp_format_vector(x$IC), "\n")
  }
  cat("Bandwidths (b):", glbfp_format_vector(x$b), "\n")
  if (!is.null(x$m)) {
    cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
  }
  if (!is.null(x$u)) {
    cat("Relative grid coordinate (u):", glbfp_format_vector(x$u), "\n")
  }
  invisible(x)
}

#' @keywords internal
glbfp_print_grid <- function(x, label) {
  cat(label, " Density Estimation on Grid:\n", sep = "")
  cat("Grid points:", nrow(x$grid), "\n")
  cat("Dimensions:", ncol(x$grid), "\n")
  cat("Grid type:", if (isTRUE(x$is_rectangular)) "rectangular" else "irregular", "\n")
  cat("Density range:", glbfp_format_vector(range(x$densities, na.rm = TRUE)), "\n")
  cat("Bandwidths (b):", glbfp_format_vector(x$b), "\n")
  if (!is.null(x$m)) {
    cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
  }
  if (!is.null(x$visited)) {
    cat("Median visited cells:", glbfp_format_number(stats::median(x$visited, na.rm = TRUE)), "\n")
  }
  if (!is.null(x$prefix_nodes)) {
    cat("Median prefix nodes:", glbfp_format_number(stats::median(x$prefix_nodes, na.rm = TRUE)), "\n")
  }
  invisible(x)
}

#' @keywords internal
glbfp_plot_grid <- function(x, contour = FALSE, type = NULL, ...) {
  if (is.null(type)) type <- if (ncol(x$grid) == 1L) "line" else if (ncol(x$grid) == 3L) "slices" else if (contour) "contour" else "surface"
  plot_density(x, type = type, ...)
}
