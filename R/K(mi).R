#' @title Function that calculates \eqn{K(m_i)}
#' @description This function computes the value of \eqn{K(m_i)}, which is an expression derived from the asymptotic properties of a statistical estimator related to the bandwidth selection in nonparametric density estimation.
#' @param mi An integer representing the index value, typically a positive integer related to the number of bins or intervals in the context of shifted histograms.
#' @return The computed value of \eqn{K(m_i)}, which is used in further statistical calculations or analyses.
#' @details The function implements the formula given by
#' \deqn{K(m_i) = \sqrt{\frac{1}{6} + \frac{1}{12m_i^2}} + \frac{4m_i^2 - 1}{6\sqrt{2}m_i} \log\left(\frac{\sqrt{3} + \sqrt{4m_i^2 + 2}}{\sqrt{4m_i^2 - 1}}\right).} This formula is derived from the theoretical analysis in the context of average shifted histograms. The terms in the formula correspond to corrections applied to the bandwidth selection to optimize statistical properties.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Calculate K for mi = 2
#'  K_mi(2)
#'
#'  # Calculate K for mi = 5
#'  K_mi(5)
#' }
#' }
#' @export
K_mi <- function(mi) {
  if (!is.numeric(mi) || length(mi) != 1L || !is.finite(mi) || mi <= 0.5) {
    stop("`mi` must be a single numeric value greater than 0.5.", call. = FALSE)
  }
  term1 <- sqrt(1 / 6 + 1 / (12 * mi^2))
  term2 <- (4 * mi^2 - 1) / (6 * sqrt(2) * mi)
  term3 <- log((sqrt(3) + sqrt(4 * mi^2 + 2)) / sqrt(4 * mi^2 - 1))

  return(term1 + term2 * term3)
}
