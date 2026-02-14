#' GLBFP density estimation on a grid
#'
#' Computes GLBFP density estimates on a regular or user-supplied grid.
#'
#' @param data Numeric matrix or data frame of observations (`n x d`).
#' @param b Positive numeric vector of bandwidths (length `d`).
#' @param m Positive integer vector of shifts (length `d`).
#' @param grid_size Integer number of grid points per dimension when
#'   `grid_points = NULL`.
#' @param grid_points Optional matrix/data frame of explicit evaluation points.
#' @param min_vals Numeric vector of lower grid bounds (length `d`).
#' @param max_vals Numeric vector of upper grid bounds (length `d`).
#'
#' @return A list with class `c("glbfp_grid", "GLBFP_estimate")` containing
#' grid coordinates, densities, uncertainty estimates, and grid metadata.
#'
#' @examples
#' b <- c(0.5, 0.5)
#' # Use a small, representative subset so examples remain fast in checks.
#' sample_data <- as.matrix(ashua[seq_len(120), -3])
#' out <- GLBFP_estimate(sample_data, b = b, m = c(1, 1), grid_size = 10)
#' out
#'
#' @export
GLBFP_estimate <- function(
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  grid_size = 20,
  grid_points = NULL,
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
) {
  data <- glbfp_validate_data(data)
  d <- ncol(data)

  b <- glbfp_validate_vector(b, d = d, name = "b", positive = TRUE)
  m <- glbfp_validate_vector(m, d = d, name = "m", positive = TRUE, integerish = TRUE)

  bounds <- glbfp_validate_bounds(min_vals, max_vals, d = d)
  min_vals <- bounds$min_vals
  max_vals <- bounds$max_vals

  grid_info <- glbfp_prepare_grid(
    data = data,
    grid_size = grid_size,
    grid_points = grid_points,
    min_vals = min_vals,
    max_vals = max_vals
  )

  results <- vapply(seq_len(nrow(grid_info$grid)), function(i) {
    res <- GLBFP(
      x = as.numeric(grid_info$grid[i, ]),
      data = data,
      b = b,
      m = m,
      min_vals = min_vals,
      max_vals = max_vals
    )
    c(estimation = res$estimation, sd = res$sd, IC_lower = res$IC[1], IC_upper = res$IC[2])
  }, numeric(4))

  result <- list(
    grid = grid_info$grid,
    densities = as.numeric(results["estimation", ]),
    sd = as.numeric(results["sd", ]),
    IC = cbind(
      IC_lower = as.numeric(results["IC_lower", ]),
      IC_upper = as.numeric(results["IC_upper", ])
    ),
    b = b,
    m = m,
    method = "GLBFP",
    grid_size = if (length(unique(grid_info$grid_dims)) == 1L) grid_info$grid_dims[1] else NA_integer_,
    grid_dims = grid_info$grid_dims,
    is_rectangular = grid_info$is_rectangular,
    col_names = grid_info$col_names
  )

  class(result) <- c("glbfp_grid", "GLBFP_estimate")
  result
}

#' @describeIn GLBFP_estimate Print method for object of class `"GLBFP_estimate"`.
#' @param x Object returned by [GLBFP_estimate()].
#' @param ... Additional arguments (unused).
#' @export
print.GLBFP_estimate <- function(x, ...) {
  cat("GLBFP Density Estimation on Grid:\n")
  cat("Grid points:", nrow(x$grid), "\n")
  cat("Dimensions:", ncol(x$grid), "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
}

#' @describeIn GLBFP_estimate Plot method for object of class `"GLBFP_estimate"`.
#' @param contour If `TRUE`, draw a contour-like 2D representation for 2D data.
#' @export
plot.GLBFP_estimate <- function(x, contour = FALSE, ...) {
  d <- ncol(x$grid)
  col_names <- x$col_names

  if (d == 1L) {
    df <- data.frame(grid = x$grid[, 1], density = x$densities)
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = grid, y = density)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "GLBFP Density Estimation", x = col_names[1], y = "Density")
    )
  }

  if (d == 2L) {
    df <- data.frame(x = x$grid[, 1], y = x$grid[, 2], z = x$densities)

    if (isTRUE(contour)) {
      if (isTRUE(x$is_rectangular)) {
        return(
          ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, z = z)) +
            ggplot2::geom_contour_filled() +
            ggplot2::labs(title = "GLBFP Density Estimation", x = col_names[1], y = col_names[2])
        )
      }
      return(
        ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = z)) +
          ggplot2::geom_point(size = 1.5) +
          ggplot2::labs(title = "GLBFP Density Estimation (Irregular Grid)", x = col_names[1], y = col_names[2], color = "Density")
      )
    }

    if (isTRUE(x$is_rectangular)) {
      x_val <- sort(unique(x$grid[, 1]))
      y_val <- sort(unique(x$grid[, 2]))
      z_val <- matrix(x$densities, nrow = length(x_val), ncol = length(y_val), byrow = FALSE)
      return(
        plotly::plot_ly(x = x_val, y = y_val, z = z_val, type = "surface") |>
          plotly::layout(scene = list(
            xaxis = list(title = col_names[1]),
            yaxis = list(title = col_names[2]),
            zaxis = list(title = "Estimated density")
          ))
      )
    }

    return(
      plotly::plot_ly(df, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers") |>
        plotly::layout(scene = list(
          xaxis = list(title = col_names[1]),
          yaxis = list(title = col_names[2]),
          zaxis = list(title = "Estimated density")
        ))
    )
  }

  stop("Plotting is only supported for dimension <= 2.", call. = FALSE)
}
