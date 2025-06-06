% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLBFP_estimate.R
\name{GLBFP_estimate}
\alias{GLBFP_estimate}
\alias{print.GLBFP_estimate}
\alias{plot.GLBFP_estimate}
\title{GLBFP Estimate for Grid}
\usage{
GLBFP_estimate(
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  grid_size = 20,
  grid_points = NULL,
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
)

\method{print}{GLBFP_estimate}(obj)

\method{plot}{GLBFP_estimate}(obj, contour = FALSE, ...)
}
\arguments{
\item{data}{Matrix or data frame of observations (n rows, d columns).}

\item{b}{Numeric vector of bandwidths for each dimension (length d).}

\item{m}{Numeric vector of shifts for each dimension (length d).}

\item{grid_size}{Integer specifying the number of grid points per dimension.}

\item{grid_points}{Numeric matrix of data framne as the grid to use to estimate LBFP estimator on each points of the grid. Default is NULL, in this case a grid is created from minimum to maximum on each dimension with step `grid_size`.}

\item{min_vals}{Minimum value of initial grid. Default value is the minimum on each column.}

\item{max_vals}{Maximum value of initial grid. Default value is the maximum on each column.}

\item{obj}{The list object returned by the `GLBFP_estimate` function.}

\item{contour}{logical (TRUE or FALSE). \code{contour = TRUE} makes a contour plot, whereas \code{contour = FALSE} makes a 3D plot using plotly.}

\item{...}{Additional arguments for customization.}
}
\value{
An object of class "GLBFP_estimate" containing the grid and density estimates.
}
\description{
Computes GLBFP density estimates on a regular grid.
}
\section{Methods (by generic)}{
\itemize{
\item \code{print(GLBFP_estimate)}: print method for object of class \code{GLBFP_estimate}
Print the density estimates calculated by the `GLBFP_estimate` function.

\item \code{plot(GLBFP_estimate)}: plot method for object of class \code{GLBFP_estimate}
Plot the GLBFP density estimate for 1D or 2D datasets.

}}
\examples{
 # Example for GLBFP_estimate function with ashua data
##  Use of grid_size
b <- c(0.5, 0.5)
out <- GLBFP_estimate(ashua[,-3], b, grid_size = 30)
out
plot(out)

## with a specific grid_points
grid_points <- data.frame(expand.grid(seq(200,250,1), seq(29,30,0.1)))
out <- GLBFP_estimate(ashua[,-3], b, grid_points = grid_points)
out
}
