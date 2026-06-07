#' LBFP density estimation on a grid
#'
#' Computes LBFP density estimates on a regular or user-supplied grid.
#'
#' @param data Numeric matrix or data frame of observations (`n x d`).
#' @param b Positive numeric vector of bandwidths (length `d`).
#' @param grid_size Integer number of grid points per dimension when
#'   `grid_points = NULL`.
#' @param grid_points Optional matrix/data frame of explicit evaluation points.
#' @param min_vals Numeric vector of lower grid bounds (length `d`).
#' @param max_vals Numeric vector of upper grid bounds (length `d`).
#'
#' @return A list with class `c("glbfp_grid", "LBFP_estimate")` containing
#' grid coordinates, densities, uncertainty estimates, and grid metadata.
#'
#' @details
#' When `grid_points` is `NULL`, a regular grid is constructed from `min_vals` to
#' `max_vals`. Custom grids may be irregular; in that case plotting uses point or
#' scatter representations instead of a surface.
#'
#' @seealso [LBFP()], [ASH_estimate()], [GLBFP_estimate()]
#'
#' @examples
#' b <- c(0.5, 0.5)
#' out <- LBFP_estimate(ashua[, -3], b = b, grid_size = 15)
#' out
#' plot(out, contour = TRUE)
#'
#' @export
LBFP_estimate <- function(
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  grid_size = 20,
  grid_points = NULL,
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
) {
  data <- glbfp_validate_data(data)
  d <- ncol(data)

  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)

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

  fast <- glbfp_evaluate_glbfp_fast(
    data = data,
    points = grid_info$grid,
    b = b,
    m = rep(1L, d),
    min_vals = min_vals,
    max_vals = max_vals
  )

  result <- list(
    grid = grid_info$grid,
    densities = as.numeric(fast$densities),
    sd = as.numeric(fast$sd),
    IC = fast$IC,
    b = b,
    method = "LBFP",
    grid_size = if (length(unique(grid_info$grid_dims)) == 1L) grid_info$grid_dims[1] else NA_integer_,
    grid_dims = grid_info$grid_dims,
    is_rectangular = grid_info$is_rectangular,
    col_names = grid_info$col_names,
    u = fast$u,
    cell_index = fast$cell_index,
    visited = fast$visited,
    prefix_nodes = fast$prefix_nodes,
    prefix_order = fast$prefix_order
  )

  class(result) <- c("glbfp_grid", "LBFP_estimate")
  result
}

#' @describeIn LBFP_estimate Print method for object of class `"LBFP_estimate"`.
#' @param x Object returned by [LBFP_estimate()].
#' @param ... Additional arguments (unused).
#' @export
print.LBFP_estimate <- function(x, ...) {
  glbfp_print_grid(x, label = "LBFP")
}

#' @describeIn LBFP_estimate Plot method for object of class `"LBFP_estimate"`.
#' @param contour If `TRUE`, draw a contour-like 2D representation for 2D data.
#' @export
plot.LBFP_estimate <- function(x, contour = FALSE, ...) {
  glbfp_plot_grid(x, contour = contour, ...)
}
