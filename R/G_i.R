#' Compute the \eqn{G(m_i)} bandwidth constant
#'
#' Computes the scalar constant \eqn{G(m_i)} used by [compute_bi_optim()].
#'
#' @param mi A single positive numeric value. In package estimators, `mi`
#'   corresponds to one component of the integer shift vector `m`.
#'
#' @return A positive numeric scalar.
#'
#' @details
#' The implemented formula is
#' \deqn{G(m_i) = \frac{1}{12}\left(1 + \frac{1}{2m_i^2}\right).}
#'
#' @seealso [compute_bi_optim()], [K_mi()], [compute_G_star()]
#'
#' @examples
#' G_i(1)
#' G_i(2)
#'
#' @export
G_i <- function(mi) {
  if (!is.numeric(mi) || length(mi) != 1L || !is.finite(mi) || mi <= 0) {
    stop("`mi` must be a single positive numeric value.", call. = FALSE)
  }
  (1 / 12) * (1 + 1 / (2 * mi^2))
}
