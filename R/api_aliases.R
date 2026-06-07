#' Lowercase aliases for the public API
#'
#' These aliases follow common R naming style while preserving the original
#' uppercase function names used in earlier versions of the package.
#'
#' @param ... Arguments passed to the corresponding uppercase function.
#'
#' @return The same object returned by the corresponding uppercase function.
#'
#' @name lowercase_aliases
NULL

#' @rdname lowercase_aliases
#' @export
ash <- function(...) {
  ASH(...)
}

#' @rdname lowercase_aliases
#' @export
lbfp <- function(...) {
  LBFP(...)
}

#' @rdname lowercase_aliases
#' @export
glbfp <- function(...) {
  GLBFP(...)
}

#' @rdname lowercase_aliases
#' @export
ash_estimate <- function(...) {
  ASH_estimate(...)
}

#' @rdname lowercase_aliases
#' @export
lbfp_estimate <- function(...) {
  LBFP_estimate(...)
}

#' @rdname lowercase_aliases
#' @export
glbfp_estimate <- function(...) {
  GLBFP_estimate(...)
}

#' @rdname lowercase_aliases
#' @export
compute_di <- function(...) {
  compute_Di(...)
}
