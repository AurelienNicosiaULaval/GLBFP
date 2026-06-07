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
#' Scott, D. W. (1992). *Multivariate Density Estimation: Theory, Practice, and
#' Visualization*. Wiley. doi:10.1002/9780470316849.
#'
#' The complete methodological citation for GLBFP has not yet been verified in
#' this repository. Add it before using this help page as publication text.
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
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
) {
  data <- glbfp_validate_data(data)
  d <- ncol(data)

  x <- glbfp_validate_vector(x, d = d, name = "x")
  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)
  m <- glbfp_validate_vector(m, d = d, name = "m", positive = TRUE, integerish = TRUE)

  bounds <- glbfp_validate_bounds(min_vals, max_vals, d = d)
  min_vals <- bounds$min_vals
  max_vals <- bounds$max_vals

  fast <- glbfp_evaluate_glbfp_fast(
    data = data,
    points = matrix(x, nrow = 1L),
    b = b,
    m = m,
    min_vals = min_vals,
    max_vals = max_vals
  )

  result <- list(
    x = x,
    estimation = fast$densities[1],
    sd = fast$sd[1],
    IC = fast$IC[1, ],
    b = b,
    m = m,
    method = "GLBFP",
    dimension = d,
    u = fast$u[1, ],
    cell_index = fast$cell_index[1, ],
    visited = fast$visited[1],
    prefix_nodes = fast$prefix_nodes[1],
    prefix_order = fast$prefix_order
  )
  class(result) <- c("glbfp_fit", "GLBFP")
  result
}

#' @describeIn GLBFP Print method for object of class `"GLBFP"`.
#' @method print GLBFP
#' @param x Object returned by [GLBFP()].
#' @param ... Additional arguments (unused).
#' @export
print.GLBFP <- function(x, ...) {
  glbfp_print_fit(x, label = "GLBFP")
}
