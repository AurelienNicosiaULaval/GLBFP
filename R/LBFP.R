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
  n <- nrow(data)

  x <- glbfp_validate_vector(x, d = d, name = "x")
  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)

  bounds <- glbfp_validate_bounds(min_vals, max_vals, d = d)
  min_vals <- bounds$min_vals
  max_vals <- bounds$max_vals

  a <- min_vals + b / 2
  max_vals <- max_vals + b / 2

  cell_count <- vapply(seq_len(d), function(i) {
    glbfp_cell_count(a[i], max_vals[i], b[i])
  }, integer(1))

  idx <- pmin(
    pmax(1L, floor((x - a) / b) + 1L),
    cell_count
  )

  lowerbound_cell <- a + (idx - 1) * b
  upperbound_cell <- a + idx * b
  mid_IK <- (lowerbound_cell + upperbound_cell) / 2

  neighbors <- as.matrix(expand.grid(rep(list(0:1), d)))

  neighbor_stats <- vapply(seq_len(nrow(neighbors)), function(j) {
    neighbor <- neighbors[j, ]
    mid <- a + (idx - 0.5) * b
    lower_bounds <- mid + (neighbor - 1) * b
    upper_bounds <- mid + neighbor * b

    count <- sum(Reduce(
      `&`,
      lapply(seq_len(d), function(i) {
        data[, i] >= lower_bounds[i] & data[, i] < upper_bounds[i]
      })
    ))

    u <- (x - (mid_IK - b / 2)) / b
    vector_c <- prod(u^neighbor * (1 - u)^(1 - neighbor))

    c(count = count, c_j = vector_c)
  }, numeric(2))

  estimation <- sum(neighbor_stats["c_j", ] * neighbor_stats["count", ]) / (n * prod(b))
  estimation <- glbfp_safe_estimation(estimation)

  sigma_hat2 <- sum(neighbor_stats["c_j", ]^2 * estimation) / (n * prod(b))
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
    method = "LBFP",
    dimension = d
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
  cat("LBFP Density Estimation:\n")
  cat("Point:", paste0("(", paste(x$x, collapse = ", "), ")"), "\n")
  cat("Estimated density:", x$estimation, "\n")
  cat("Estimated standard error:", x$sd, "\n")
  cat("95% confidence interval:", paste(x$IC, collapse = ", "), "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
}
