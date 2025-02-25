#' GLBFP Density Estimator at a Single Point
#'
#' This function computes the GLBFP (General Linear Blend Frequency Polygon) estimate of the density
#' at a single given point. The GLBFP is a continuous density estimator inspired by the frequency polygon.
#'
#' @param x Numeric vector representing the point where the density is to be estimated.
#' @param data Numeric matrix or data frame representing the dataset. Each row is an observation.
#' @param b Numeric vector representing the bin width for each dimension. Default is optimal bi calculated using \code{compute_bi_optim}.
#' @param m Numeric vector representing the shift for each dimension. Default is 1 on each dimension (similar to LBFP estimation).
#' @param min_vals Minimum value of initial grid (I_0(k)). Default value is the minimum on each column.
#' @param max_vals Maximum value of initial grid (I_0(k)). Default value is the maximum on each column.
#' @examples
#' # Example for GLBFP function with bivariate normale distribution
#' data <- as.data.frame(mvrnorm(n = 1000, mu =  c(0, 0), Sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)))
#' x <- c(0.5, 0.5)
#' ## with a selected bin width
#' b <- c(0.5, 0.5)
#' m <- c(1,1)
#' GLBFP(x, data, b, m)
#' ## with optimal b and shift m=3
#' GLBFP(x, data, m = c(3,3))
#' @return A list containing the estimated density at point `x` using the GLBFP estimator, the bin width vector `b`, the shift vector `m`, the point of estimation `x`, the standard deviation of the estimation and finally 95% confidence interval.
#' @import data.table
#' @export
GLBFP <- function(x, data, b = compute_bi_optim(data, m = rep(1,  ncol(data))), m = rep(1,  ncol(data)), min_vals =  apply(data, 2, min), max_vals =  apply(data, 2, max) ) {
  # Input validation
  if (!is.numeric(x) || !is.numeric(b) || !is.numeric(m) || (!is.matrix(data) && !is.data.frame(data))) {
    stop("Inputs must be numeric vectors for x, m, and b, and a matrix or data frame for data.")
  }
  d <- ncol(data)
  if (length(x) != d || length(b) != d || length(m) != d || length(min_vals) != d || length(max_vals) != d) {
    stop("The point, bin width, min_vals, max_val and m vectors must match the dataset dimensions.")
  }
  if (any(b <= 0)) stop("All bin widths must be positive.")
  if (any(m <= 0)) stop("All m values must be positive.")
  
  # Convert data to matrix for faster processing
  n <- nrow(data)
  data <- as.matrix(data)
  
  # Compute the size of subcells (delta) within each dimension
  delta <- b / m
  
  # Compute minimum and maximum bounds for each dimension
  a <- min_vals + delta/2
  max_vals <- max_vals + delta/2
  
  # Compute the index of the cell containing the point x
  idx <- pmin(pmax(1, floor((x - a) / delta) + 1), sapply(1:d, function(i) length(seq(a[i], max_vals[i], by = delta[i])) - 1))
  
  
  # Compute the lower and upper bounds of the current cell
  lowerbound_cell <- a + (idx - 1) * delta
  upperbound_cell <- a + idx * delta
  
  
  
  # Generate all possible neighbors (binary combinations)
  neighbors <- expand.grid(rep(list(0:1), d))
  
  # Generate all possible local indices for m
  local_indices_matrix <- expand.grid(lapply(m, function(mi) seq(1 - mi, mi - 1)))
  local_indices_matrix <- as.matrix(local_indices_matrix)
  
  # Compute the midpoint of the current cell
  mid <- a + (idx - 0.5) * delta
  u <- (x - (mid-delta/2)) / delta
  
  # Pre-compute weights for all combinations of local indices
  weights <- apply(local_indices_matrix, 1, function(row) {
    prod(pmax(0, 1 - abs(row) / m))
  })
  
  # Compute ash estimation for each neighbor
  counts <- apply(neighbors, 1, function(neighbor) {
    # Compute the starting point of the subcell for the current neighbor
    x0 <- lowerbound_cell + neighbor * delta
    
    # Compute bounds for each dimension and all local indices
    cell_bounds <- lapply(1:d, function(j) {
      lower <- x0[j] + (local_indices_matrix[, j]  - 0.5) * delta[j]
      upper <- x0[j] + (local_indices_matrix[, j]  + 0.5) * delta[j]
      list(lower = lower, upper = upper)
    })
    
    # Filter valid points within the bounds for each dimension
    valid_points <- lapply(1:d, function(j) {
      outer(data[, j], cell_bounds[[j]]$lower, `>=`) &
        outer(data[, j], cell_bounds[[j]]$upper, `<`)
    })
    
    # Combine results across dimensions and count valid rows
    valid_combined <- Reduce(`&`, valid_points)
    counts <- colSums(valid_combined)
    
    # Compute ash estimation
    ash_estimation <- sum(weights * counts) / (n * prod(b))
    
    # Compute c_j, the contribution of this neighbor
    vector_c <- prod(u^neighbor * (1 - u)^(1 - neighbor))
    return(c(ash_estimation, vector_c))
  })
  
  # Final estimation by combining ash estimations and weights
  estimation <- sum(counts[2, ] * counts[1, ])
  
  # ensure estimation is greater than 0 (boundary)
  if(estimation<0 || is.nan(estimation) ) estimation <- 0
  
  
  # variance estimation
  sigma_hat2 <- 1 / (n * prod(b)) * prod((2*m^2+1-6*u*(1-u))/(3*m^2))*estimation
  
  
  # confiance interval
  IC = c(estimation + qnorm(0.025)/sqrt((n * prod(b)))*sqrt(sigma_hat2), estimation + qnorm(0.975)/sqrt((n * prod(b)))*sqrt(sigma_hat2))
  
  # Return the result as a list
  result <- list(
    x = x,
    estimation = estimation,
    sd = sqrt(sigma_hat2),
    IC = IC,
    b = b,
    m = m
  )
  class(result) <- "GLBFP"
  return(result)
}


#' @describeIn GLBFP print method for object of class "GLBFP"
#' @method print GLBFP
#' @export
print.GLBFP <- function(obj, ...) {
  cat("GLBFP Density Estimation:\n")
  cat("Point of estimation:", paste0("(", paste(obj$x, collapse = ", "), ")"), "\n")
  cat("Estimated density:", obj$estimation, "\n")
  cat("Estimated standard error",obj$sd, "\n")
  cat("95% confidence interval:",obj$IC,"\n")
  cat("Bandwidths (b):", paste(obj$b, collapse = ", "), "\n")
  cat("Shifts (m):", paste(obj$m, collapse = ", "), "\n")
}
