#' ASH density estimation on a grid
#'
#' Computes ASH density estimates on a regular or user-supplied grid.
#'
#' @param data Numeric matrix or data frame of observations (`n x d`).
#' @param b Positive numeric vector of bandwidths (length `d`).
#' @param m Positive integer vector of shifts (length `d`).
#' @param grid_size Integer number of grid points per dimension when
#'   `grid_points = NULL`.
#' @param grid_points Optional matrix/data frame of explicit evaluation points.
#' @param min_vals Numeric vector of lower grid bounds (length `d`).
#' @param max_vals Numeric vector of upper grid bounds (length `d`).
#'
#' @return A list with class `c("glbfp_grid", "ASH_estimate")` containing
#' grid coordinates, densities, and grid metadata.
#'
#' @details
#' When `grid_points` is `NULL`, a regular grid is constructed from `min_vals` to
#' `max_vals`. Custom grids may be irregular; in that case plotting uses point or
#' scatter representations instead of a surface.
#'
#' @seealso [ASH()], [LBFP_estimate()], [GLBFP_estimate()]
#'
#' @examples
#' b <- c(0.5, 0.5)
#' # Use a small, representative subset so examples remain fast in checks.
#' sample_data <- as.matrix(ashua[seq_len(120), -3])
#' out <- ASH_estimate(sample_data, b = b, m = c(1, 1), grid_size = 10)
#' out
#'
#' @export
ASH_estimate <- function(
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  grid_size = 20,
  grid_points = NULL,
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
) {
  data <- glbfp_validate_data(data)
  d <- ncol(data)

  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)
  m <- glbfp_validate_vector(m, d = d, name = "m", positive = TRUE, integerish = TRUE)

  bounds <- glbfp_validate_bounds(min_vals, max_vals, d = d)
  min_vals <- bounds$min_vals
  max_vals <- bounds$max_vals

  grid_info <- glbfp_prepare_grid(
    data = data,
    grid_size = grid_size,
    grid_points = grid_points,
    min_vals = min_vals,
    max_vals = max_vals
  )

  fast <- glbfp_evaluate_ash_fast(
    data = data,
    points = grid_info$grid,
    b = b,
    m = m,
    min_vals = min_vals,
    max_vals = max_vals
  )

  result <- list(
    grid = grid_info$grid,
    densities = as.numeric(fast$densities),
    b = b,
    m = m,
    method = "ASH",
    grid_size = if (length(unique(grid_info$grid_dims)) == 1L) grid_info$grid_dims[1] else NA_integer_,
    grid_dims = grid_info$grid_dims,
    is_rectangular = grid_info$is_rectangular,
    col_names = grid_info$col_names,
    cell_index = fast$cell_index,
    visited = fast$visited,
    prefix_nodes = fast$prefix_nodes,
    prefix_order = fast$prefix_order
  )

  class(result) <- c("glbfp_grid", "ASH_estimate")
  result
}

#' @describeIn ASH_estimate Print method for object of class `"ASH_estimate"`.
#' @param x Object from [ASH_estimate()] to print.
#' @param ... Additional arguments (unused).
#' @method print ASH_estimate
#' @export
print.ASH_estimate <- function(x, ...) {
  glbfp_print_grid(x, label = "ASH")
}

#' @describeIn ASH_estimate Plot method for object of class `"ASH_estimate"`.
#' @param contour If `TRUE`, draw a contour-like 2D representation for 2D data.
#' @method plot ASH_estimate
#' @export
plot.ASH_estimate <- function(x, contour = FALSE, ...) {
  glbfp_plot_grid(x, contour = contour, ...)
}
