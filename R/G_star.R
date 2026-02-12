#' @title Function that calculates \eqn{G^*}
#' @description This function computes the value of \eqn{G^*}, a constant derived from the asymptotic analysis of bandwidth selection in high-dimensional settings.
#' @param d An integer representing the dimension of the dataset.
#' @return The computed value of \eqn{G^*}, used in the context of high-dimensional nonparametric density estimation.
#' @details The function implements the formula given by
#' \eqn{G^* = 2^{\frac{3(d - 4)}{2(4 + d)}} \cdot \exp\left(\frac{1}{4 + d}\right) \cdot \frac{\pi^{d/2}}{4 + d}}. This formula is derived from the analysis of high-dimensional density estimation techniques, where it acts as a scaling factor in the calculation of bandwidth.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Calculate G* for d = 2
#'  compute_G_star(2)
#'
#'  # Calculate G* for d = 5
#'  compute_G_star(5)
#' }
#' }
#' @export
compute_G_star <- function(d) {
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d < 1 || abs(d - round(d)) > sqrt(.Machine$double.eps)) {
    stop("`d` must be a positive integer.", call. = FALSE)
  }
  d <- as.integer(round(d))
  two_term <- 2^((3 * (d - 4)) / (2 * (4 + d)))
  e1_term <- exp(1 / (4 + d))
  pi_term <- pi^(d / 2) / (4 + d)
  return(two_term * e1_term * pi_term)
}
