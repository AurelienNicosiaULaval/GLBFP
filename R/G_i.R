#' @title Function that calculates \eqn{G(m_i)}
#' @description This function computes the value of \eqn{G(m_i)}, which is a term related to the asymptotic properties of the variance of the average shifted histogram estimator.
#' @param mi An integer representing the index value, typically a positive integer related to the number of bins or intervals in the context of shifted histograms.
#' @return The computed value of \eqn{G(m_i)}, used in the context of variance calculations in nonparametric density estimation.
#' @details The function implements the formula given by
#' \eqn{G(m_i) = \frac{1}{12} \left(1 + \frac{1}{2m_i^2}\right)}. This expression is derived from theoretical considerations on the bias and variance trade-off in density estimation using average shifted histograms.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Calculate G for mi = 2
#'  G_i(2)
#'
#'  # Calculate G for mi = 5
#'  G_i(5)
#' }
#' }
#' @export
G_i <- function(mi) {
  if (!is.numeric(mi) || length(mi) != 1L || !is.finite(mi) || mi <= 0) {
    stop("`mi` must be a single positive numeric value.", call. = FALSE)
  }
  return((1 / 12) * (1 + 1 / (2 * mi^2)))
}
