#' Summarize GLBFP fit objects
#'
#' Summarizes objects returned by [ASH()], [LBFP()], [GLBFP()] and
#' their grid counterparts.
#'
#' @param object A fitted object of class `"glbfp_fit"` or `"glbfp_grid"`.
#' @param ... Additional arguments (unused).
#'
#' @return A list with class `"summary.glbfp_fit"` or `"summary.glbfp_grid"`.
#' @export
summary.glbfp_fit <- function(object, ...) {
  out <- list(
    method = object$method,
    dimension = object$dimension,
    point = object$x,
    estimation = object$estimation,
    sd = object$sd %||% NA_real_,
    IC = object$IC %||% c(NA_real_, NA_real_)
  )
  class(out) <- "summary.glbfp_fit"
  out
}

#' @export
print.summary.glbfp_fit <- function(x, ...) {
  cat("Method:", x$method, "\n")
  cat("Dimension:", x$dimension, "\n")
  cat("Point:", paste(x$point, collapse = ", "), "\n")
  cat("Estimation:", x$estimation, "\n")
  if (!is.na(x$sd)) {
    cat("Standard error:", x$sd, "\n")
  }
  if (all(is.finite(x$IC))) {
    cat("95% CI:", paste(x$IC, collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
summary.glbfp_grid <- function(object, ...) {
  out <- list(
    method = object$method,
    dimension = ncol(object$grid),
    n_grid = nrow(object$grid),
    density_min = min(object$densities),
    density_median = stats::median(object$densities),
    density_mean = mean(object$densities),
    density_max = max(object$densities)
  )
  class(out) <- "summary.glbfp_grid"
  out
}

#' @export
print.summary.glbfp_grid <- function(x, ...) {
  cat("Method:", x$method, "\n")
  cat("Dimension:", x$dimension, "\n")
  cat("Grid points:", x$n_grid, "\n")
  cat("Density range:", paste(c(x$density_min, x$density_max), collapse = " to "), "\n")
  cat("Density median:", x$density_median, "\n")
  cat("Density mean:", x$density_mean, "\n")
  invisible(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
