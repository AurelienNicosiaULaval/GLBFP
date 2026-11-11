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
#' @param engine Computational engine; see [prepare_glbfp()].
#' @param inference Default `"none"` returns `NA` uncertainty fields. `"iid"`
#'   requires explicit fixed bandwidths and anchor and independent observations.
#'   The sample variance of observation weights estimates the variance of the
#'   estimate. Normal 95 percent intervals target the expected smoothed density,
#'   without bias correction or guaranteed finite-sample coverage. Evaluation
#'   points must also be fixed independently of the observations.
#'   `"stationary_mixing"` implements the LBFP plug-in asymptotic interval of
#'   Carbon and Duchesne (2026), requiring a complete rectangular indexed
#'   field, explicit fixed bandwidths and anchor, and fixed query points. It
#'   applies to GLBFP only when all `m` equal 1. Stationarity, mixing, regularity,
#'   asymptotic bandwidth conditions and negligible bias are caller assumptions;
#'   the function cannot verify them. The interval targets the marginal density,
#'   is pointwise and asymptotic, and does not guarantee finite-sample coverage.
#'   A warning is issued: strong dependence can cause severe undercoverage in
#'   finite samples. Zero plug-in density gives an unavailable interval.
#' @section Grid conventions:
#' Fine counting cells are half-open with width `b/m` and lower-bound anchor
#' `min_vals`. Interpolation is between their centers. Unoccupied cells have
#' zero counts on an unbounded lattice. `max_vals` controls the display grid,
#' not the density support. See [prepare_glbfp()] for exact conventions.
#' `data` may also be a [glbfp_field()]; index coordinates are kept separate.
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
  b = NULL,
  grid_size = 20,
  grid_points = NULL,
  min_vals = NULL,
  max_vals = NULL,
  engine = c("auto", "R", "cpp_direct", "cpp_prefix"),
  inference = c("none", "iid", "stationary_mixing")
) {
  glbfp_fit_grid(
    data, grid_size, grid_points, b, NULL, "LBFP",
    min_vals, max_vals, match.arg(engine), match.arg(inference)
  )
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
