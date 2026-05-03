#' Compute the \eqn{G^*} bandwidth constant
#'
#' Computes the dimension-dependent constant \eqn{G^*} used by
#' [compute_bi_optim()].
#'
#' @param d A single positive integer giving the data dimension.
#'
#' @return A positive numeric scalar.
#'
#' @details
#' The implemented formula is
#' \deqn{G^* =
#' 2^{\frac{3(d - 4)}{2(4 + d)}}
#' \exp\left(\frac{1}{4 + d}\right)
#' \frac{\pi^{d / 2}}{4 + d}.}
#'
#' @seealso [compute_bi_optim()], [G_i()], [K_mi()]
#'
#' @examples
#' compute_G_star(1)
#' compute_G_star(2)
#'
#' @export
compute_G_star <- function(d) {
  if (
    !is.numeric(d) ||
      length(d) != 1L ||
      !is.finite(d) ||
      d < 1 ||
      abs(d - round(d)) > sqrt(.Machine$double.eps)
  ) {
    stop("`d` must be a positive integer.", call. = FALSE)
  }
  d <- as.integer(round(d))
  two_term <- 2^((3 * (d - 4)) / (2 * (4 + d)))
  e1_term <- exp(1 / (4 + d))
  pi_term <- pi^(d / 2) / (4 + d)
  two_term * e1_term * pi_term
}
