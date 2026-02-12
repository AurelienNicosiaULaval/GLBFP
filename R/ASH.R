#' Averaged Shifted Histogram (ASH) estimator at a single point
#'
#' Computes the ASH density estimate at point `x`.
#'
#' @param x Numeric vector with one coordinate per data dimension.
#' @param data Numeric matrix or data frame of observations (`n x d`).
#' @param b Positive numeric vector of bandwidths (length `d`).
#' @param m Positive integer vector of shifts (length `d`).
#' @param min_vals Numeric vector of lower grid bounds (length `d`).
#' @param max_vals Numeric vector of upper grid bounds (length `d`).
#'
#' @return A list with class `c("glbfp_fit", "ASH")` containing:
#' `x`, `estimation`, `b`, `m`, `method`, and `dimension`.
#'
#' @examples
#' x <- c(200, 30)
#' b <- c(0.5, 0.5)
#' m <- c(1, 1)
#' ASH(x, ashua[, -3], b = b, m = m)
#'
#' @export
ASH <- function(
  x,
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
) {
  data <- glbfp_validate_data(data)
  d <- ncol(data)
  n <- nrow(data)

  x <- glbfp_validate_vector(x, d = d, name = "x")
  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)
  m <- glbfp_validate_vector(m, d = d, name = "m", positive = TRUE, integerish = TRUE)

  bounds <- glbfp_validate_bounds(min_vals, max_vals, d = d)
  min_vals <- bounds$min_vals
  max_vals <- bounds$max_vals

  delta <- b / m
  a <- min_vals

  cell_count <- vapply(seq_len(d), function(i) {
    glbfp_cell_count(a[i], max_vals[i], delta[i])
  }, integer(1))

  idx <- pmin(
    pmax(1L, floor((x - a) / delta) + 1L),
    cell_count
  )

  lowerbound_cell <- a + (idx - 1) * delta
  upperbound_cell <- a + idx * delta
  x0 <- (lowerbound_cell + upperbound_cell) / 2

  local_indices_matrix <- as.matrix(expand.grid(lapply(m, function(mi) seq(1 - mi, mi - 1))))

  weights <- apply(local_indices_matrix, 1, function(l) {
    prod(pmax(0, 1 - abs(l) / m))
  })

  ash_estimation <- sum(vapply(seq_len(nrow(local_indices_matrix)), function(index) {
    l <- local_indices_matrix[index, ]

    valid_points <- Reduce(
      `&`,
      lapply(seq_len(d), function(j) {
        lower <- x0[j] + (l[j] - 0.5) * delta[j]
        upper <- x0[j] + (l[j] + 0.5) * delta[j]
        data[, j] >= lower & data[, j] < upper
      })
    )

    weights[index] * sum(valid_points)
  }, numeric(1)))

  estimation <- glbfp_safe_estimation(ash_estimation / (n * prod(b)))

  result <- list(
    x = x,
    estimation = estimation,
    b = b,
    m = m,
    method = "ASH",
    dimension = d
  )
  class(result) <- c("glbfp_fit", "ASH")
  result
}

#' @describeIn ASH Print method for object of class `"ASH"`.
#' @method print ASH
#' @param x Object of class `"ASH"`.
#' @param ... Additional arguments (unused).
#' @export
print.ASH <- function(x, ...) {
  cat("ASH Density Estimation:\n")
  cat("Point:", paste0("(", paste(x$x, collapse = ", "), ")"), "\n")
  cat("Estimated density:", x$estimation, "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
}
