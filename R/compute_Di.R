#' Leave-one-out self-support scores for grid density estimators
#'
#' Computes the fixed-grid leave-one-out score
#' \deqn{D_i = 1 - \widehat f_{(-i)}(X_i) / \widehat f(X_i)}
#' for observations in `data`. The grid, bandwidths, and estimator parameters
#' are held fixed when the contribution of observation `i` is removed.
#'
#' @param data Numeric matrix or data frame of observations (`n x d`).
#' @param b Positive numeric vector of bandwidths (length `d`).
#' @param m Positive integer vector of shifts (length `d`). Ignored for
#'   `estimator = "LBFP"`, where `m = 1`.
#' @param estimator Character string. One of `"GLBFP"`, `"LBFP"`, or `"ASH"`.
#' @param min_vals Numeric vector of lower grid bounds (length `d`).
#' @param max_vals Numeric vector of upper grid bounds (length `d`).
#'
#' @return A list with class `"glbfp_di"` containing the score vector `D`,
#' fitted densities, leave-one-out densities, self-weights, and metadata.
#'
#' @examples
#' x <- as.matrix(ashua[seq_len(80), -3])
#' b <- c(0.5, 0.5)
#' out <- compute_Di(x, b = b, m = c(1, 1), estimator = "GLBFP")
#' summary(out)
#'
#' @export
compute_Di <- function(
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  estimator = c("GLBFP", "LBFP", "ASH"),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
) {
  data <- glbfp_validate_data(data)
  d <- ncol(data)
  n <- nrow(data)
  if (n < 2L) {
    stop("`data` must contain at least two observations.", call. = FALSE)
  }

  estimator <- match.arg(estimator)
  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)
  if (identical(estimator, "LBFP")) {
    m <- rep(1L, d)
  } else {
    m <- glbfp_validate_vector(m, d = d, name = "m", positive = TRUE, integerish = TRUE)
  }

  bounds <- glbfp_validate_bounds(min_vals, max_vals, d = d)
  min_vals <- bounds$min_vals
  max_vals <- bounds$max_vals

  delta <- b / m
  self_state <- glbfp_count_state(data = data, delta = delta, min_vals = min_vals)

  if (identical(estimator, "ASH")) {
    fast <- glbfp_evaluate_ash_fast(
      data = data,
      points = data,
      b = b,
      m = m,
      min_vals = min_vals,
      max_vals = max_vals,
      self_keys = self_state$data_keys
    )
  } else {
    fast <- glbfp_evaluate_glbfp_fast(
      data = data,
      points = data,
      b = b,
      m = m,
      min_vals = min_vals,
      max_vals = max_vals,
      self_keys = self_state$data_keys
    )
  }

  volume <- prod(b)
  D <- glbfp_compute_di_from_density(
    density = fast$densities,
    self_weight = fast$self_weight,
    n = n,
    volume = volume
  )
  density_loo <- (n * fast$densities * volume - fast$self_weight) / ((n - 1) * volume)
  density_loo[!is.finite(density_loo)] <- NA_real_

  result <- list(
    D = D,
    D_positive = pmax(D, 0),
    density = fast$densities,
    density_loo = density_loo,
    self_weight = fast$self_weight,
    b = b,
    m = m,
    method = estimator,
    dimension = d,
    n = n,
    min_vals = min_vals,
    max_vals = max_vals,
    cell_index = fast$cell_index,
    visited = fast$visited,
    prefix_nodes = fast$prefix_nodes,
    prefix_order = fast$prefix_order
  )
  class(result) <- "glbfp_di"
  result
}

#' @export
print.glbfp_di <- function(x, ...) {
  cat("Leave-one-out D_i scores\n")
  cat("Method:", x$method, "\n")
  cat("Observations:", x$n, "\n")
  cat("Dimension:", x$dimension, "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  if (!identical(x$method, "LBFP")) {
    cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
  }
  cat("D_i range:", paste(range(x$D, na.rm = TRUE), collapse = " to "), "\n")
  invisible(x)
}

#' @export
summary.glbfp_di <- function(object, ...) {
  probs <- c(0, 0.25, 0.5, 0.75, 1)
  out <- list(
    method = object$method,
    n = object$n,
    dimension = object$dimension,
    D = stats::quantile(object$D, probs = probs, na.rm = TRUE),
    D_mean = mean(object$D, na.rm = TRUE),
    D_missing = sum(is.na(object$D)),
    density_min = min(object$density, na.rm = TRUE),
    density_median = stats::median(object$density, na.rm = TRUE),
    density_max = max(object$density, na.rm = TRUE),
    visited_median = stats::median(object$visited, na.rm = TRUE),
    prefix_nodes_median = stats::median(object$prefix_nodes, na.rm = TRUE)
  )
  class(out) <- "summary.glbfp_di"
  out
}

#' @export
print.summary.glbfp_di <- function(x, ...) {
  cat("D_i score summary\n")
  cat("Method:", x$method, "\n")
  cat("Observations:", x$n, "\n")
  cat("Dimension:", x$dimension, "\n")
  cat("D_i quantiles:\n")
  print(x$D)
  cat("D_i mean:", x$D_mean, "\n")
  cat("Missing D_i:", x$D_missing, "\n")
  cat("Density range:", paste(c(x$density_min, x$density_max), collapse = " to "), "\n")
  cat("Median visited cells:", x$visited_median, "\n")
  cat("Median prefix nodes:", x$prefix_nodes_median, "\n")
  invisible(x)
}

#' @export
plot.glbfp_di <- function(x, type = c("index", "density"), positive = FALSE, ...) {
  type <- match.arg(type)
  df <- as.data.frame(x)
  df$score <- if (isTRUE(positive)) df$D_positive else df$D

  if (identical(type, "index")) {
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = observation, y = score)) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        ggplot2::geom_point(size = 1.4, alpha = 0.8) +
        ggplot2::labs(
          title = paste(x$method, "leave-one-out D_i scores"),
          x = "Observation",
          y = if (isTRUE(positive)) "Positive part of D_i" else "D_i"
        ) +
        ggplot2::theme_minimal()
    )
  }

  ggplot2::ggplot(df, ggplot2::aes(x = density, y = score)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
    ggplot2::geom_point(size = 1.4, alpha = 0.8) +
    ggplot2::labs(
      title = paste(x$method, "D_i scores by fitted density"),
      x = "Fitted density",
      y = if (isTRUE(positive)) "Positive part of D_i" else "D_i"
    ) +
    ggplot2::theme_minimal()
}
