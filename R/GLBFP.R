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
  n <- nrow(data)

  x <- glbfp_validate_vector(x, d = d, name = "x")
  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)
  m <- glbfp_validate_vector(m, d = d, name = "m", positive = TRUE, integerish = TRUE)

  bounds <- glbfp_validate_bounds(min_vals, max_vals, d = d)
  min_vals <- bounds$min_vals
  max_vals <- bounds$max_vals

  delta <- b / m
  a <- min_vals + delta / 2
  max_vals <- max_vals + delta / 2

  cell_count <- vapply(seq_len(d), function(i) {
    glbfp_cell_count(a[i], max_vals[i], delta[i])
  }, integer(1))

  idx <- pmin(
    pmax(1L, floor((x - a) / delta) + 1L),
    cell_count
  )

  lowerbound_cell <- a + (idx - 1) * delta

  neighbors <- as.matrix(expand.grid(rep(list(0:1), d)))
  local_indices_matrix <- as.matrix(expand.grid(lapply(m, function(mi) seq(1 - mi, mi - 1))))

  mid <- a + (idx - 0.5) * delta
  u <- (x - (mid - delta / 2)) / delta

  weights <- apply(local_indices_matrix, 1, function(row) {
    prod(pmax(0, 1 - abs(row) / m))
  })

  counts <- vapply(seq_len(nrow(neighbors)), function(i) {
    neighbor <- neighbors[i, ]
    x0 <- lowerbound_cell + neighbor * delta

    valid_combined <- Reduce(`&`, lapply(seq_len(d), function(j) {
      lower <- x0[j] + (local_indices_matrix[, j] - 0.5) * delta[j]
      upper <- x0[j] + (local_indices_matrix[, j] + 0.5) * delta[j]
      outer(data[, j], lower, `>=`) & outer(data[, j], upper, `<`)
    }))

    cell_counts <- colSums(valid_combined)
    ash_estimation <- sum(weights * cell_counts) / (n * prod(b))
    vector_c <- prod(u^neighbor * (1 - u)^(1 - neighbor))

    c(ash = ash_estimation, c_j = vector_c)
  }, numeric(2))

  estimation <- sum(counts["c_j", ] * counts["ash", ])
  estimation <- glbfp_safe_estimation(estimation)

  sigma_hat2 <- (1 / (n * prod(b))) *
    prod((2 * m^2 + 1 - 6 * u * (1 - u)) / (3 * m^2)) *
    estimation
  sigma_hat2 <- glbfp_safe_estimation(sigma_hat2)

  se_term <- sqrt(sigma_hat2) / sqrt(n * prod(b))
  IC <- c(
    estimation + stats::qnorm(0.025) * se_term,
    estimation + stats::qnorm(0.975) * se_term
  )

  result <- list(
    x = x,
    estimation = estimation,
    sd = sqrt(sigma_hat2),
    IC = IC,
    b = b,
    m = m,
    method = "GLBFP",
    dimension = d
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
  cat("GLBFP Density Estimation:\n")
  cat("Point:", paste0("(", paste(x$x, collapse = ", "), ")"), "\n")
  cat("Estimated density:", x$estimation, "\n")
  cat("Estimated standard error:", x$sd, "\n")
  cat("95% confidence interval:", paste(x$IC, collapse = ", "), "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
}
