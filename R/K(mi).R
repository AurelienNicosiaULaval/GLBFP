#' Compute the \eqn{K(m_i)} bandwidth constant
#'
#' Computes the scalar constant \eqn{K(m_i)} used by [compute_bi_optim()].
#'
#' @param mi A single numeric value greater than 0.5. In package estimators,
#'   `mi` corresponds to one component of the integer shift vector `m`.
#'
#' @return A positive numeric scalar.
#'
#' @details
#' The implemented formula is
#' \deqn{
#' K(m_i) =
#' \sqrt{\frac{1}{6} + \frac{1}{12m_i^2}} +
#' \frac{4m_i^2 - 1}{6\sqrt{2}m_i}
#' \log\left(
#' \frac{\sqrt{3} + \sqrt{4m_i^2 + 2}}{\sqrt{4m_i^2 - 1}}
#' \right).
#' }
#'
#' @seealso [compute_bi_optim()], [G_i()], [compute_G_star()]
#'
#' @examples
#' K_mi(1)
#' K_mi(2)
#'
#' @export
K_mi <- function(mi) {
  if (!is.numeric(mi) || length(mi) != 1L || !is.finite(mi) || mi <= 0.5) {
    stop("`mi` must be a single numeric value greater than 0.5.", call. = FALSE)
  }
  inverse <- 1 / mi
  term1 <- sqrt(1 / 6 + inverse^2 / 12)
  arg <- sqrt(3) * inverse / sqrt(4 - inverse^2)
  ratio <- if (arg == 0) 1 else asinh(arg) / arg
  term2 <- sqrt(12 - 3 * inverse^2) / (6 * sqrt(2)) * ratio
  term1 + term2
}
