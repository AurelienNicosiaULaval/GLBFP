#' @title Function that computes the optimal \eqn{b_i}
#' @description This function calculates the optimal values of \eqn{b_i}, which are used in nonparametric density estimation in multiple dimensions. The calculation involves the covariance matrix of the input data, the inverse covariance matrix, and various other terms derived from the dimensionality of the data and statistical properties of the estimator.
#' @param data A data frame representing the input dataset, where rows correspond to observations and columns to variables (or dimensions).
#' @param m A vector of integers representing the values of \eqn{m_i} for each dimension of the dataset. Default value is 1 on each dimension.
#' @return A vector containing the optimal \eqn{b_i} values for each dimension.
#' @details This function performs the following steps:
#' - Computes the covariance matrix \eqn{\Sigma} and its inverse and determinant.
#' - Calculates \eqn{G_0} as the product of \eqn{K(m_i)} for each dimension.
#' - Computes \eqn{G_i} values for each dimension.
#' - Calculates \eqn{G^*} using the formula that depends on the dimensionality \eqn{d}.
#' - Uses these intermediate results to compute the optimal \eqn{b_i} values, following the formula:
#' \deqn{b_i = G^* \cdot G_0^{2 / (4 + d)} \cdot \prod(G_i)^{1 / (2 \cdot (4 + d))} \cdot (\det(\Sigma))^{1 / (2 \cdot (4 + d))} \cdot (\Sigma^{-1}_{ii})^{1 / (2 \cdot (4 + d))} \cdot G_i^{-1 / 2} \cdot (\Sigma^{-1}_{ii})^{-1 / (2 \cdot (4 + d))} \cdot n^{-1 / (4 + d)}}.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Example usage with a sample dataset
#'  sample_data <- data.frame(x = rnorm(100), y = rnorm(100))
#'  m_values <- c(2, 2)
#'  compute_bi_optim(sample_data, m_values)
#' }
#' }
#' @export
compute_bi_optim <- function(data, m = rep(1, ncol(data))) {

  # Check if the input is a data.frame
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data should be a data.frame.")
  }

  n <- nrow(data)
  d <- ncol(data)

  # Estimate the variance-covariance matrix Σ and its inverse
  Sigma <- cov(data)
  Sigma_inv <- solve(Sigma)
  Sigma_det <- det(Sigma)  # Calculate the determinant of Σ

  # Calculate G_0 as the product of K(m_i)
  G_0 <- prod(sapply(m, K_mi))

  G_i_values <- sapply(m, G_i)  # G_i is calculated for each m_i

  G_star <- compute_G_star(d)

  # Vectorized calculation of b_i for each dimension
  bi <- G_star * G_0^(2 / (4 + d)) * prod(G_i_values)^(1 / (2 * (4 + d))) *
    (Sigma_det^(1 / (2 * (4 + d)))) *
    (prod(diag(Sigma_inv))^(1 / (2 * (4 + d)))) *
    G_i_values^(-1 / 2) *
    diag(Sigma_inv)^(-1 / (2 * (4 + d))) *
    n^(-1 / (4 + d))

  return(bi)
}
