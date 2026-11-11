#' General Linear Blend Frequency Polygon (GLBFP) estimator at a single point
#'
#' Computes the GLBFP density estimate at point `x`.
#'
#' @param x Numeric vector with one coordinate per data dimension.
#' @param data Numeric matrix or data frame of observations (`n x d`).
#' @param b Positive numeric vector of bandwidths (length `d`).
#' @param m Positive integer vector of shifts (length `d`).
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
#' @return A list with class `c("glbfp_fit", "GLBFP")` containing:
#' `x`, `estimation`, `sd`, `IC`, `b`, `m`, `method`, and `dimension`.
#'
#' @details
#' `GLBFP()` generalizes the linear blend frequency polygon workflow through the
#' positive integer shift vector `m`. Missing and non-finite values are not
#' accepted; remove or impute them before calling the estimator.
#'
#' @seealso [GLBFP_estimate()], [ASH()], [LBFP()], [compute_bi_optim()]
#'
#' @references
#' Carbon, M., and Duchesne, T. (2026; online 2025). Asymptotic normality of
#' multivariate frequency polygons for stationary random fields. *Annals of the
#' Institute of Statistical Mathematics*, 78, 297-326.
#' doi:10.1007/s10463-025-00952-x.
#'
#' Scott, D. W. (1992). *Multivariate Density Estimation: Theory, Practice, and
#' Visualization*. Wiley. doi:10.1002/9780470316849.
#'
#' The precise GLBFP construction follows equations (2)-(7) of Carbon (2024),
#' internal working note, Le fameux GLBFP. Its Gaussian-reference bandwidth
#' rule is heuristic; asymptotic normality for this generalization is not
#' established by that note.
#'
#' @examples
#' x <- c(200, 30)
#' b <- c(0.5, 0.5)
#' m <- c(1, 1)
#' GLBFP(x, ashua[, -3], b = b, m = m)
#'
#' @export
GLBFP <- function(
  x,
  data,
  b = NULL,
  m = NULL,
  min_vals = NULL,
  max_vals = NULL,
  engine = c("auto", "R", "cpp_direct", "cpp_prefix"),
  inference = c("none", "iid", "stationary_mixing")
) {
  glbfp_fit_point(
    x, data, b, m, "GLBFP",
    min_vals, max_vals, match.arg(engine), match.arg(inference)
  )
}

#' @describeIn GLBFP Print method for object of class `"GLBFP"`.
#' @method print GLBFP
#' @param x Object returned by [GLBFP()].
#' @param ... Additional arguments (unused).
#' @export
print.GLBFP <- function(x, ...) {
  glbfp_print_fit(x, label = "GLBFP")
}
