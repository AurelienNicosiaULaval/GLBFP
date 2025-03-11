#' LBFP Density Estimator at a Single Point
#'
#' This function computes the LBFP (Linear Blend Frequency Polygon) estimate of the density
#' at a single given point. The LBFP is a continuous density estimator inspired by the frequency polygon.
#'
#' @param x Numeric vector representing the point where the density is to be estimated.
#' @param data Numeric matrix or data frame representing the dataset. Each row is an observation.
#' @param b Numeric vector representing the bin width for each dimension. Default is optimal bi calculated using \code{compute_bi_optim}.
#' @param min_vals Minimum value of initial grid (I_0(k)). Default value is the minimum on each column.
#' @param max_vals Maximum value of initial grid (I_0(k)). Default value is the maximum on each column.
#' @examples
#' # Example for GLBFP function with ashua data
#' x <- c(200, 30)
#' ## with a selected bin width
#' b <- c(0.5, 0.5)
#' LBFP(x, ashua[,-3], b)
#' ## with optimal b 
#' LBFP(x, ashua[,-3])
#' @return A list containing the estimated density at point x using the LBFP estimator, the bin width vector b, the standard deviation of the estimation and finally 95% confidence interval.
#' @importFrom data.table as.data.table
#' @export
LBFP <- function(x, data, b = compute_bi_optim(data, m = rep(1,  ncol(data))), min_vals =  apply(data, 2, min), max_vals =  apply(data, 2, max) ) {
  # Ensure input is valid
  if (!is.numeric(x) || !is.numeric(b) || (!is.matrix(data) && !is.data.frame(data))) {
    stop("Inputs must be numeric vectors for x and b, and a matrix or data frame for data.")
  }
  d <- ncol(data)
  if (length(x) != d || length(b) != d || length(min_vals) != d || length(max_vals) != d) {
    stop("The point, min_vals, max_vals and bin width vectors must match the dataset dimensions.")
  }
  if (any(b <= 0)) {
    stop("All bin widths must be positive.")
  }

  n <- nrow(data)

  # Convert to data.table for faster operations
  data <- data.table::as.data.table(data)

  # Calculate bounds for the hypercube containing x
  a <- min_vals + b/2 # Lower bounds
  max_vals <- max_vals+ b/2 # Upper bounds

  idx <- sapply(1:d, function(i) {
    index <- floor((x[i] - a[i]) / b[i]) + 1
    # Ensure idx stays within valid limits
    index <- pmin(pmax(1, index), length(seq(a[i], max_vals[i], by = b[i])) - 1)
    return(index)
  })

 # Compute of I(k)
 lowerbound_cell_IK <- sapply(1:d, function(i) { a[i] + (idx[i] - 1) * b[i] })
 upperbound_cell_IK <- sapply(1:d, function(i) { a[i] + idx[i] * b[i] })
 mid_IK <- (lowerbound_cell_IK + upperbound_cell_IK)/2

  # Calculate the number of points and c_j in each neighboring hypercube
  neighbors <- expand.grid(rep(list(0:1), d))

  counts <- sapply(1:nrow(neighbors), function(j) {
    neighbor <- as.numeric(neighbors[j, ])
    mid <- a + (idx-1/2)*b
    lower_bounds <- mid + (neighbor-1) * b
    upper_bounds <- mid+ neighbor* b

    # Use data.table for fast filtering
    count <- data[ , sum(Reduce(`&`, lapply(1:d, function(i) {
      data[[i]] >= lower_bounds[i] & data[[i]] < upper_bounds[i]
    })))]

    # compute c_j
    u <- (x - (mid_IK - b/2))/b
    vector_c <-prod(u^(neighbor)*(1-u)^(1-neighbor))
    out = rbind(count, vector_c)
    return(out)
  })

  # Finalize the estimation
  estimation <- sum(counts[2,] * counts[1,]) / (n * prod(b))

  # ensure estimation is greater than 0 (boundary)
  if(estimation<0 || is.nan(estimation) ) estimation <- 0

  # variance estimation
  sigma_hat2 <- sum(counts[2,]^2*estimation)/(n*prod(b))

  # confiance interval
  IC = c(estimation + qnorm(0.025)/sqrt((n * prod(b)))*sqrt(sigma_hat2), estimation + qnorm(0.975)/sqrt((n * prod(b)))*sqrt(sigma_hat2))

  result <- list(
    x = x,
    estimation = estimation,
    sd = sqrt(sigma_hat2) ,
    IC = IC,
    b = b
  )
  class(result) <- "LBFP"
  return(result)
}


#' @describeIn LBFP Print object of class \code{"LBFP"}
#' @method print LBFP
#' @param obj Object form \code{LBFP} to print
#' @export
print.LBFP <- function(obj) {
  cat("LBFP Density Estimation:\n")
  cat("Point of estimation:", paste0("(", paste(obj$x, collapse = ", "), ")"), "\n")
  cat("Estimated density:", obj$estimation,"\n")
  cat("Estimated standard error",obj$sd, "\n")
  cat("95% confidence interval:",obj$IC,"\n")
  cat("Bandwidths (b):", paste(obj$b, collapse = ", "), "\n")
}
