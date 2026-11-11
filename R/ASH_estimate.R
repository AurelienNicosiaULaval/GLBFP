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
#' @param engine Computational engine; see [prepare_glbfp()].
#' @param inference Default `"none"` returns `NA` uncertainty fields. `"iid"`
#'   requires explicit fixed bandwidths and anchor and independent observations.
#'   The sample variance of observation weights estimates the variance of the
#'   estimate. Normal 95 percent intervals target the expected smoothed density,
#'   without bias correction or guaranteed finite-sample coverage. Evaluation
#'   points must also be fixed independently of the observations. No dependent
#'   field interval is implemented for ASH.
#' @section Grid conventions:
#' Fine counting cells are half-open with width `b/m` and lower-bound anchor
#' `min_vals`. Interpolation is between their centers. Unoccupied cells have
#' zero counts on an unbounded lattice. `max_vals` controls the display grid,
#' not the density support. See [prepare_glbfp()] for exact conventions.
#' `data` may also be a [glbfp_field()]; index coordinates are kept separate.
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
  b = NULL,
  m = NULL,
  grid_size = 20,
  grid_points = NULL,
  min_vals = NULL,
  max_vals = NULL,
  engine = c("auto", "R", "cpp_direct", "cpp_prefix"),
  inference = c("none", "iid")
) {
  glbfp_fit_grid(
    data, grid_size, grid_points, b, m, "ASH",
    min_vals, max_vals, match.arg(engine), match.arg(inference)
  )
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
