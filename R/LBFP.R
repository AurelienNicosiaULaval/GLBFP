#' Linear Blend Frequency Polygon (LBFP) estimator at a single point
#'
#' Computes the LBFP density estimate at point `x`.
#'
#' @param x Numeric vector with one coordinate per data dimension.
#' @param data Numeric matrix or data frame of observations (`n x d`).
#' @param b Positive numeric vector of bandwidths (length `d`).
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
#'   field interval is implemented.
#' @section Grid conventions:
#' Fine counting cells are half-open with width `b/m` and lower-bound anchor
#' `min_vals`. Interpolation is between their centers. Unoccupied cells have
#' zero counts on an unbounded lattice. `max_vals` controls the display grid,
#' not the density support. See [prepare_glbfp()] for exact conventions.
#' `data` may also be a [glbfp_field()]; index coordinates are kept separate.
#' @return A list with class `c("glbfp_fit", "LBFP")` containing:
#' `x`, `estimation`, `sd`, `IC`, `b`, `method`, and `dimension`.
#'
#' @details
#' The estimate is obtained by linear blending of neighboring histogram bin
#' heights. Missing and non-finite values are not accepted; remove or impute them
#' before calling the estimator.
#'
#' @seealso [LBFP_estimate()], [ASH()], [GLBFP()], [compute_bi_optim()]
#'
#' @references
#' Scott, D. W. (1992). *Multivariate Density Estimation: Theory, Practice, and
#' Visualization*. Wiley. doi:10.1002/9780470316849.
#'
#' Terrell, G. R., and Scott, D. W. (1985). Oversmoothed Nonparametric Density
#' Estimates. *Journal of the American Statistical Association*, 80(389),
#' 209-214. doi:10.1080/01621459.1985.10477163.
#'
#' @examples
#' x <- c(200, 30)
#' b <- c(0.5, 0.5)
#' LBFP(x, ashua[, -3], b = b)
#'
#' @export
LBFP <- function(
  x,
  data,
  b = NULL,
  min_vals = NULL,
  max_vals = NULL,
  engine = c("auto", "R", "cpp_direct", "cpp_prefix"),
  inference = c("none", "iid")
) {
  glbfp_fit_point(
    x, data, b, NULL, "LBFP",
    min_vals, max_vals, match.arg(engine), match.arg(inference)
  )
}

#' @describeIn LBFP Print method for object of class `"LBFP"`.
#' @method print LBFP
#' @param x Object from [LBFP()].
#' @param ... Additional arguments (unused).
#' @export
print.LBFP <- function(x, ...) {
  glbfp_print_fit(x, label = "LBFP")
}
