#' LBFP Density Estimator
#'
#' This function computes the LBFP density estimate for every observation in the dataset.
#'
#' @param data Numeric matrix or data frame representing the dataset. Each row is an observation.
#' @param b Numeric vector representing the bin width for each dimension. Default is optimal bi calculated using \code{compute_bi_optim}.
#' @param grid_size Size of the grid for each dimension. Default is 20.
#' @param grid_points Numeric matrix of data framne as the grid to use to estimate LBFP estimator on each points of the grid. Default is NULL, in this case a grid is created from minimum to maximum on each dimension with step `grid_size`.
#' @param min_vals Minimum value of initial grid. Default value is the minimum on each column.
#' @param max_vals Maximum value of initial grid. Default value is the maximum on each column.
#' @examples
#'  # Example for LBFP_estimate function with ashua data
#' ##  Use of grid_size
#' b <- c(0.5, 0.5)
#' out <- LBFP_estimate(ashua[,-3], b, grid_size = 30)
#' out
#' plot(out)
#'
#' ## with a specific grid_points
#' grid_points <- data.frame(expand.grid(seq(200,250,1), seq(29,30,0.1)))
#' out <- LBFP_estimate(ashua[,-3], b, grid_points = grid_points)
#' out
#' @importFrom data.table as.data.table
#' @return A list containing the estimated density values for each point in the grid, the bin width vector `b`, the grid points and the `grid_size`.
#' @export
LBFP_estimate <- function(data, b = compute_bi_optim(data, m = rep(1,  ncol(data))), grid_size = 20, grid_points = NULL, min_vals =  apply(data, 2, min), max_vals =  apply(data, 2, max) ) {
  # Ensure input is valid
  d <- ncol(data)
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame.")
  }
  if (!is.numeric(b) || length(b) != ncol(data)|| length(min_vals) != d || length(max_vals) != d ) {
    stop("Bin width, min_vals, max_vals and m vector must be numeric and match the number of columns in data.")
  }
  if (any(b <= 0)) {
    stop("All bin widths must be positive.")
  }
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
    res <- LBFP(point, data, b, min_vals, max_vals)
    c(estimation = res$estimation, sd = res$sd, IC_lower = res$IC[1], IC_upper = res$IC[2])
  })
  
  results_mat <- do.call(rbind, results)
  
  result <- list(grid = grid_points,
                 densities = results_mat[, "estimation"],
                 sd = results_mat[, "sd"],
                 IC = results_mat[, c("IC_lower", "IC_upper")],
                 b = b,
                 grid_size = grid_size)
  
  class(result) <- "LBFP_estimate"
  return(result)
}


#' @describeIn LBFP_estimate Print object of class \code{"LBFP_estimate"}
#'
#' Print the density estimates calculated by the `LBFP_estimate` function.
#'
#' @param x The list object returned by the `LBFP_estimate` function.
#' @param ... Additional arguments (unused).
#' @export
print.LBFP_estimate <- function(x, ...) {
  cat("LBFP Density Estimation on Grid:\n")
  cat("Grid size:", x$grid_size, "\n")
  cat("Bandwidths (b):", paste(x$b, collapse = ", "), "\n")
  cat("Number of estimated densities:", length(x$densities), "\n")
}

#' @describeIn LBFP_estimate Plot object of class \code{"LBFP_estimate"}
#'
#' Plot the LBFP density estimate for 1D or 2D datasets.
#'
#' @param x The list object returned by the `LBFP_estimate` function.
#' @param contour logical (TRUE or FALSE). \code{contour = TRUE} makes a contour plot, whereas \code{contour = FALSE} makes a 3D plot using plotly.
#' @param ... Additional arguments for customization.
#' @importFrom magrittr %>%
#' @export
plot.LBFP_estimate <- function(x, contour = FALSE, ... ) {
  if (!is.list(x) || !"densities" %in% names(x) || !"grid" %in% names(x)) {
    stop("Invalid LBFP_estimate object. Expected a list with 'estimations' and 'grid' components.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting. Please install it.")
  }
  data <- x$grid
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("grid must be a matrix or data frame.")
  }
  d <- ncol(data)
  if (d > 2) {
    stop("Plotting is only available for 1D or 2D data.")
  }

  if (d == 1) {
    # 1D plot of density estimates using ggplot2
    df <- data.frame(Values = x$grid, Density = x$densities)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Values, y = Density)) +
      ggplot2::geom_line(color = "blue") +
      ggplot2::geom_point(color = "red") +
      ggplot2::geom_area(fill = "lightblue", alpha = 0.5) +
      ggplot2::labs(title = "Linear Blend Frequency Polygon (LBFP) en 1D",
                    x = colnames(data)[1],
                    y = "Estimated density") +
      ggplot2::theme_minimal()
    print(p)
  } else if (d == 2) {
    # 2D contour or 3D plot using ggplot2 or plotly
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("The plotly package is required for 3D plotting. Please install it.")
    }
    if (contour) {
      df <- data.frame(x = x$grid[, 1], y = x$grid[, 2], z = x$densities)
      # Contour plot using ggplot2
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, z = z)) +
        ggplot2::geom_contour_filled() +
        ggplot2::labs(title = "LBFP Contour Plot (2D)",
                      x = colnames(data)[1],
                      y = colnames(data)[2],
                      fill = "Density estimation") +
        ggplot2::theme_minimal()
      print(p)
    } else {
      # 3D Plot using plotly
      x_val <- unique(x$grid[,1])
      y_val <- unique(x$grid[,2])
      z_val <-matrix(x$densities, nrow = x$grid_size, ncol = x$grid_size)
      plotly::plot_ly() %>%
        plotly::add_surface(x = ~x_val, y = ~y_val, z = ~z_val) %>%
        plotly::layout(scene = list(
                 xaxis = list(title = colnames(data)[1]),
                 yaxis = list(title = colnames(data)[2]),
                 zaxis = list(title = "Estimated density")
               ))
       }
  }
}
