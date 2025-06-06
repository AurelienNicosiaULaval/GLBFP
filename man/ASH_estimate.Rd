% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ASH_estimate.R
\name{ASH_estimate}
\alias{ASH_estimate}
\alias{print.ASH_estimate}
\alias{plot.ASH_estimate}
\title{ASH Estimate for Grid}
\usage{
ASH_estimate(
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  grid_size = 20,
  grid_points = NULL,
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
)

\method{print}{ASH_estimate}(obj, ...)

\method{plot}{ASH_estimate}(obj, contour = FALSE, ...)

\method{print}{ASH_estimate}(obj, ...)

\method{plot}{ASH_estimate}(obj, contour = FALSE, ...)
}
\arguments{
\item{data}{Matrix or data frame of observations (n rows, d columns).}

\item{b}{Numeric vector of bandwidths for each dimension (length d).}

\item{m}{Numeric vector of shifts for each dimension (length d).}

\item{grid_size}{Integer specifying the number of grid points per dimension.}

\item{grid_points}{Numeric matrix of data framne as the grid to use to estimate LBFP estimator on each points of the grid. Default is NULL, in this case a grid is created from minimum to maximum on each dimension with step `grid_size`.}

\item{min_vals}{Minimum value of initial grid. Default value is the minimum on each column.}

\item{max_vals}{Maximum value of initial grid. Default value is the maximum on each column.}

\item{obj}{Object form \code{ASH_estimate} to plot}

\item{...}{Additional arguments for customization.}

\item{contour}{Logical to plot contour plot or interactive 3D plot. Default is FALSE (3D plot)}
}
\value{
An object of class "ASH_estimate" containing the grid and density estimates.
}
\description{
Computes ASH density estimates on a regular grid.
}
\section{Methods (by generic)}{
\itemize{
\item \code{print(ASH_estimate)}: Print method for object of class \code{"ASH_estimate"}
Print the density estimates calculated by the `ASH_estimate` function.

\item \code{plot(ASH_estimate)}: Plot method for object of class \code{"ASH_estimate"}
Plot the ASH density estimate for 1D or 2D datasets.

\item \code{print(ASH_estimate)}: Print object of class \code{"ASH_estimate"}

\item \code{plot(ASH_estimate)}: Plot object of class \code{"ASH_estimate"}

}}
\examples{
 # Example for ASH_estimate function with ashua data
##  Use of grid_size
b <- c(0.5, 0.5)
out <- ASH_estimate(ashua[,-3], b, grid_size = 30)
out
plot(out)

## with a specific grid_points
grid_points <- data.frame(expand.grid(seq(200,250,1), seq(29,30,0.1)))
out <- ASH_estimate(ashua[,-3], b, grid_points = grid_points)
out
}
