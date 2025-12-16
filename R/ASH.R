#' Averaged Shifted Histogram (ASH) Estimator
#'
#' Computes the ASH density estimate at a single point `x`.
#'
#' @param x Numeric vector representing the point at which to estimate density.
#' @param data Matrix or data frame of observations (n rows, d columns).
#' @param b Numeric vector representing the bin width for each dimension. Default is optimal bi calculated using \code{compute_bi_optim}.
#' @param m Numeric vector representing the shift for each dimension. Default is 1 on each dimension.
#' @param min_vals Minimum value of initial grid. Default value is the minimum on each column.
#' @param max_vals Maximum value of initial grid. Default value is the maximum on each column.
#' @return A list containing the estimated density at point `x` using the ASH estimator, the bin width vector `b`, the shift vector `m` and `x`
#' @examples
#' # Example for ASH function with ashua data
#' x <- c(200, 30)
#' ## with a selected bin width
#' b <- c(0.5, 0.5)
#' m <- c(1,1)
#' ASH(x, ashua[,-3], b, m)
#' ## with optimal b and shift m=3
#' ASH(x, ashua[,-3], m = c(3,3))
#' @export
#'
ASH <- function(
  x,
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
) {
  # Input validation
  if (
    !is.numeric(x) ||
      !is.numeric(b) ||
      !is.numeric(m) ||
      (!is.matrix(data) && !is.data.frame(data))
  ) {
    stop(
      "Inputs must be numeric vectors for x, m, and b, and a matrix or data frame for data."
    )
  }
  d <- ncol(data)
  if (
    length(x) != d ||
      length(b) != d ||
      length(m) != d ||
      length(min_vals) != d ||
      length(max_vals) != d
  ) {
    stop(
      "The point, bin width, min_vals, max_vals and m vectors must match the dataset dimensions."
    )
  }
  if (any(b <= 0)) {
    stop("All bin widths must be positive.")
  }
  if (any(m <= 0)) {
    stop("All m values must be positive.")
  }

  # Convert data to matrix for faster processing
  n <- nrow(data)
  data <- as.matrix(data)

  # Compute the minimum, maximum, and index for each dimension
  a <- min_vals #+ b/2
  max_vals <- max_vals #+ b/2

  # compute delta
  delta <- b / m

  idx <- pmin(
    pmax(1, floor((x - a) / delta) + 1),
    sapply(1:d, function(i) length(seq(a[i], max_vals[i], by = delta[i])) - 1)
  )

  # Compute cell bounds and the center of the current cell
  lowerbound_cell <- a + (idx - 1) * delta
  upperbound_cell <- a + idx * delta
  x0 <- (lowerbound_cell + upperbound_cell) / 2

  # generate all local indices

  local_indices_matrix <- as.matrix(expand.grid(lapply(m, function(mi) {
    seq(1 - mi, mi - 1)
  })))

  # Precompute weights for all local indices
  weights <- apply(local_indices_matrix, 1, function(l) {
    prod(pmax(0, 1 - abs(l) / m))
  })

  # Compute ash estimation using vectorized operations
  ash_estimation <- sum(sapply(1:nrow(local_indices_matrix), function(index) {
    l <- local_indices_matrix[index, ]
    # Compute bounds for each dimension
    cell_bounds <- lapply(1:d, function(j) {
      lower <- x0[j] + (l[j] - 0.5) * delta[j]
      upper <- x0[j] + (l[j] + 0.5) * delta[j]
      list(lower = lower, upper = upper)
    })

    # Vectorized filtering for valid points
    valid_points <- Reduce(
      `&`,
      lapply(1:d, function(j) {
        outer(data[, j], cell_bounds[[j]]$lower, `>=`) &
          outer(data[, j], cell_bounds[[j]]$upper, `<`)
      })
    )
    count_in_cell <- colSums(valid_points)

    # Compute weighted contribution for this cell
    weights[index] * sum(count_in_cell)
  }))

  # Normalize ash estimation
  ash_estimation <- ash_estimation / (n * prod(b))
  # ensure estimation is greater than 0 (boundary)
  if (ash_estimation < 0 || is.nan(ash_estimation)) {
    ash_estimation <- 0
  }

  # Return the result as a list
  result <- list(
    x = x,
    estimation = ash_estimation,
    b = b,
    m = m
  )
  class(result) <- "ASH"
  return(result)
}

#' @describeIn ASH Print method for object of class \code{"ASH"}
#' @method print ASH
#' @param x Object of class \code{"ASH"}.
#' @param ... Additional arguments (unused).
#' @export
#' @keywords internal
print.ASH <- function(x, ...) {
  cat("ASH Density Estimation:\n")
  cat(
    "Point of estimation:",
    paste0("(", paste(x$x, collapse = ", "), ")"),
    "\n"
  )
  cat("Estimated density:", x$estimation, "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
}
