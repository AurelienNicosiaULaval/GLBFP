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
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
) {
  data <- glbfp_validate_data(data)
  d <- ncol(data)

  x <- glbfp_validate_vector(x, d = d, name = "x")
  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)

  bounds <- glbfp_validate_bounds(min_vals, max_vals, d = d)
  min_vals <- bounds$min_vals
  max_vals <- bounds$max_vals

  fast <- glbfp_evaluate_glbfp_fast(
    data = data,
    points = matrix(x, nrow = 1L),
    b = b,
    m = rep(1L, d),
    min_vals = min_vals,
    max_vals = max_vals
  )

  result <- list(
    x = x,
    estimation = fast$densities[1],
    sd = fast$sd[1],
    IC = fast$IC[1, ],
    b = b,
    method = "LBFP",
    dimension = d,
    u = fast$u[1, ],
    cell_index = fast$cell_index[1, ],
    visited = fast$visited[1],
    prefix_nodes = fast$prefix_nodes[1],
    prefix_order = fast$prefix_order
  )
  class(result) <- c("glbfp_fit", "LBFP")
  result
}

#' @describeIn LBFP Print method for object of class `"LBFP"`.
#' @method print LBFP
#' @param x Object from [LBFP()].
#' @param ... Additional arguments (unused).
#' @export
print.LBFP <- function(x, ...) {
  glbfp_print_fit(x, label = "LBFP")
}
