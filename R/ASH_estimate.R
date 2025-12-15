

#' @title ASH Estimate for Grid
#' @description Computes ASH density estimates on a regular grid.
#' @param data Matrix or data frame of observations (n rows, d columns).
#' @param b Numeric vector of bandwidths for each dimension (length d).
#' @param m Numeric vector of shifts for each dimension (length d).
#' @param grid_size Integer specifying the number of grid points per dimension.
#' @param grid_points Numeric matrix of data framne as the grid to use to estimate LBFP estimator on each points of the grid. Default is NULL, in this case a grid is created from minimum to maximum on each dimension with step `grid_size`.
#' @param min_vals Minimum value of initial grid. Default value is the minimum on each column.
#' @param max_vals Maximum value of initial grid. Default value is the maximum on each column.
#' @examples
#'  # Example for ASH_estimate function with ashua data
#' ##  Use of grid_size
#' b <- c(0.5, 0.5)
#' out <- ASH_estimate(ashua[,-3], b, grid_size = 30)
#' out
#' plot(out)
#'
#' ## with a specific grid_points
#' grid_points <- data.frame(expand.grid(seq(200,250,1), seq(29,30,0.1)))
#' out <- ASH_estimate(ashua[,-3], b, grid_points = grid_points)
#' out
#' @return An object of class "ASH_estimate" containing the grid and density estimates.
#' @export
ASH_estimate <- function(data, b = compute_bi_optim(data, m = rep(1,  ncol(data))), m = rep(1,  ncol(data)), grid_size = 20, grid_points = NULL, min_vals =  apply(data, 2, min), max_vals =  apply(data, 2, max) ) {
  # Input validation
  if ( !is.numeric(b) || !is.numeric(m) || (!is.matrix(data) && !is.data.frame(data))) {
    stop("Inputs must be numeric vectors for x, m, and b, and a matrix or data frame for data.")
  }
  d <- ncol(data)
  if ( length(b) != d || length(m) != d || length(min_vals) != d || length(max_vals) != d ) {
    stop("The point, bin width, min_vals, max_vals and m vectors must match the dataset dimensions.")
  }
  if (any(b <= 0)) stop("All bin widths must be positive.")
  if (any(m <= 0)) stop("All m values must be positive.")

  if ((grid_size <= 0)) {
    stop("Grid size must be positive.")
  }
  if (!is.null(grid_points) && (!is.matrix(grid_points) && !is.data.frame(grid_points) || ncol(grid_points)!= ncol(data))){
    stop("Grid_points must be numeric and match the number of columns in the data")
  }
  if(!is.null(grid_points))   grid_size <- length(unique(grid_points[,1]))

  if(is.null(grid_points)){
    # Get the bounds for the grid
    d <- ncol(data)
    grid_points <- vector("list", d)
    # Convert to data.table for faster operations
    data <- data.table::as.data.table(data)
    # Calculate the bounds for each dimension
    for (i in 1:d) {
      grid_points[[i]] <- seq(min(data[[i]]), max(data[[i]]), length.out = grid_size)
    }

    # Create the grid
    grid_points <- as.matrix(expand.grid(grid_points))
  }

  # Compute estimation, sd, and IC on each grid point
  results <- lapply(1:nrow(grid_points), function(i) {
    point <- as.numeric(grid_points[i, ])
    res <- ASH(point, data, b, m, min_vals, max_vals)
    c(estimation = res$estimation)
  })
  
  results_mat <- do.call(rbind, results)
  
  result <- list(grid = grid_points,
                 densities = results_mat[, "estimation"],
                 b = b,
                 m = m,
                 grid_size = grid_size)
  
  class(result) <- "ASH_estimate"
  return(result)
}

#' @describeIn ASH_estimate Print object of class \code{"ASH_estimate"}
#' @param x Object form \code{ASH_estimate} to print
#' @param ... Additional arguments (unused).
#' @method print ASH_estimate
#' @export
print.ASH_estimate <- function(x, ...) {
  cat("ASH Density Estimation on Grid:\n")
  cat("Grid size:", x$grid_size, "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  cat("Shifts (m):", paste(x$m, collapse = ", "), "\n")
  cat("Number of estimated densities:", length(x$densities), "\n")
}

#' @describeIn ASH_estimate Plot object of class \code{"ASH_estimate"}
#' @param x Object form \code{ASH_estimate} to plot
#' @param contour Logical to plot contour plot or interactive 3D plot. Default is FALSE (3D plot)
#' @method plot ASH_estimate
#' @export
plot.ASH_estimate <- function(x, contour = FALSE, ...) {
  d <- ncol(x$grid)
  if (d == 1) {
    df <- data.frame(grid = x$grid[, 1], density = x$densities)
    ggplot2::ggplot(df,  ggplot2::aes(x = grid, y = density)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = "ASH Density Estimation", x = "Grid", y = "Density")
  } else if (d == 2) {
    df <- data.frame(x = x$grid[, 1], y = x$grid[, 2], z = x$densities)
    if (contour) {
      ggplot2::ggplot(df,  ggplot2::aes(x = x, y = y, z = z)) +
        ggplot2::geom_contour_filled() +
        ggplot2::labs(title = "ASH Density Estimation",  x = colnames(data)[1],
                      y = colnames(data)[2])
    } else {
      x_val <- unique(x$grid[,1])
      y_val <- unique(x$grid[,2])
      z_val <-  matrix(x$densities, nrow = x$grid_size, ncol = x$grid_size)
      plotly::plot_ly() %>%
        plotly::add_surface(x = ~x_val, y = ~y_val, z = ~z_val) %>%
        plotly::layout(scene = list(
          xaxis = list(title = colnames(data)[1]),
          yaxis = list(title = colnames(data)[2]),
          zaxis = list(title = "Estimated density")
        ))
    }
  } else {
    stop("Plotting is only supported for d <= 2.")
  }
}
