#' Compute bandwidth vector \eqn{b_i}
#'
#' Computes a plug-in bandwidth vector used by GLBFP/LBFP/ASH estimators.
#' The function validates numeric inputs, stabilizes near-singular covariance
#' matrices with a small ridge if needed, and returns strictly positive
#' bandwidths.
#'
#' @param data A numeric matrix or data frame where rows are observations and
#'   columns are variables.
#' @param m A positive integer vector of shifts, one value per dimension.
#'
#' @return A numeric vector of positive bandwidths with one value per column in
#'   `data`.
#'
#' @examples
#' set.seed(1)
#' x <- cbind(rnorm(200), rnorm(200))
#' compute_bi_optim(x, m = c(1, 1))
#'
#' @export
compute_bi_optim <- function(data, m = rep(1, ncol(data))) {
  data <- glbfp_validate_data(data)
  d <- ncol(data)
  n <- nrow(data)

  m <- glbfp_validate_vector(m, d = d, name = "m", positive = TRUE, integerish = TRUE)

  Sigma <- stats::cov(data)
  Sigma_stable <- glbfp_stabilize_covariance(Sigma)
  Sigma_inv <- Sigma_stable$Sigma_inv
  Sigma_det <- Sigma_stable$Sigma_det

  G_0 <- prod(vapply(m, K_mi, numeric(1)))
  G_i_values <- vapply(m, G_i, numeric(1))
  G_star <- compute_G_star(d)

  bi <- G_star * G_0^(2 / (4 + d)) * prod(G_i_values)^(1 / (2 * (4 + d))) *
    Sigma_det^(1 / (2 * (4 + d))) *
    prod(diag(Sigma_inv))^(1 / (2 * (4 + d))) *
    G_i_values^(-1 / 2) *
    diag(Sigma_inv)^(-1 / (2 * (4 + d))) *
    n^(-1 / (4 + d))

  if (any(!is.finite(bi)) || any(bi <= 0)) {
    stop("Failed to compute strictly positive finite bandwidths.", call. = FALSE)
  }

  unname(as.numeric(bi))
}
