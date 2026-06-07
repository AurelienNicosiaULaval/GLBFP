# Shared print and plot helpers.

#' @keywords internal
glbfp_format_number <- function(x, digits = 6L) {
  format(signif(x, digits), trim = TRUE)
}

#' @keywords internal
glbfp_format_vector <- function(x, digits = 6L) {
  paste(glbfp_format_number(x, digits = digits), collapse = ", ")
}

#' @keywords internal
glbfp_print_fit <- function(x, label) {
  cat(label, " Density Estimation:\n", sep = "")
  cat("Point:", paste0("(", glbfp_format_vector(x$x), ")"), "\n")
  cat("Estimated density:", glbfp_format_number(x$estimation), "\n")
  if (!is.null(x$sd) && is.finite(x$sd)) {
    cat("Estimated standard error:", glbfp_format_number(x$sd), "\n")
  }
  if (!is.null(x$IC) && all(is.finite(x$IC))) {
    cat("95% confidence interval:", glbfp_format_vector(x$IC), "\n")
  }
  cat("Bandwidths (b):", glbfp_format_vector(x$b), "\n")
  if (!is.null(x$m)) {
    cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
  }
  if (!is.null(x$u)) {
    cat("Relative grid coordinate (u):", glbfp_format_vector(x$u), "\n")
  }
  invisible(x)
}

#' @keywords internal
glbfp_print_grid <- function(x, label) {
  cat(label, " Density Estimation on Grid:\n", sep = "")
  cat("Grid points:", nrow(x$grid), "\n")
  cat("Dimensions:", ncol(x$grid), "\n")
  cat("Grid type:", if (isTRUE(x$is_rectangular)) "rectangular" else "irregular", "\n")
  cat("Density range:", glbfp_format_vector(range(x$densities, na.rm = TRUE)), "\n")
  cat("Bandwidths (b):", glbfp_format_vector(x$b), "\n")
  if (!is.null(x$m)) {
    cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
  }
  if (!is.null(x$visited)) {
    cat("Median visited cells:", glbfp_format_number(stats::median(x$visited, na.rm = TRUE)), "\n")
  }
  if (!is.null(x$prefix_nodes)) {
    cat("Median prefix nodes:", glbfp_format_number(stats::median(x$prefix_nodes, na.rm = TRUE)), "\n")
  }
  invisible(x)
}

#' @keywords internal
glbfp_plot_grid <- function(x, contour = FALSE, ...) {
  d <- ncol(x$grid)
  col_names <- x$col_names
  title <- paste(x$method, "Density Estimation")

  if (d == 1L) {
    df <- data.frame(grid = x$grid[, 1], density = x$densities)
    df <- df[order(df$grid), , drop = FALSE]
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = grid, y = density)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = title, x = col_names[1], y = "Estimated density") +
        ggplot2::theme_minimal()
    )
  }

  if (d == 2L) {
    df <- data.frame(x = x$grid[, 1], y = x$grid[, 2], z = x$densities)

    if (isTRUE(contour)) {
      if (isTRUE(x$is_rectangular)) {
        return(
          ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, z = z)) +
            ggplot2::geom_contour_filled() +
            ggplot2::labs(title = title, x = col_names[1], y = col_names[2], fill = "Density") +
            ggplot2::theme_minimal()
        )
      }
      return(
        ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = z)) +
          ggplot2::geom_point(size = 1.5) +
          ggplot2::labs(
            title = paste(title, "(Irregular Grid)"),
            x = col_names[1],
            y = col_names[2],
            color = "Density"
          ) +
          ggplot2::theme_minimal()
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
