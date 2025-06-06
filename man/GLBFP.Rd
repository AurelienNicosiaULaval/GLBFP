% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLBFP.R
\name{GLBFP}
\alias{GLBFP}
\alias{print.GLBFP}
\title{GLBFP Density Estimator at a Single Point}
\usage{
GLBFP(
  x,
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
)

\method{print}{GLBFP}(obj, ...)
}
\arguments{
\item{x}{Numeric vector representing the point where the density is to be estimated.}

\item{data}{Numeric matrix or data frame representing the dataset. Each row is an observation.}

\item{b}{Numeric vector representing the bin width for each dimension. Default is optimal bi calculated using \code{compute_bi_optim}.}

\item{m}{Numeric vector representing the shift for each dimension. Default is 1 on each dimension (similar to LBFP estimation).}

\item{min_vals}{Minimum value of initial grid (I_0(k)). Default value is the minimum on each column.}

\item{max_vals}{Maximum value of initial grid (I_0(k)). Default value is the maximum on each column.}
}
\value{
A list containing the estimated density at point `x` using the GLBFP estimator, the bin width vector `b`, the shift vector `m`, the point of estimation `x`, the standard deviation of the estimation and finally 95% confidence interval.
}
\description{
This function computes the GLBFP (General Linear Blend Frequency Polygon) estimate of the density
at a single given point. The GLBFP is a continuous density estimator inspired by the frequency polygon.
}
\section{Methods (by generic)}{
\itemize{
\item \code{print(GLBFP)}: print method for object of class "GLBFP"

}}
\examples{
# Example for GLBFP function with ashua data
x <- c(200, 30)
## with a selected bin width
b <- c(0.5, 0.5)
m <- c(1,1)
GLBFP(x, ashua[,-3], b, m)
## with optimal b and shift m=3
GLBFP(x, ashua[,-3], m = c(3,3))
}
