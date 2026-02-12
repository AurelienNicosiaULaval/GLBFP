#' Predict from GLBFP fit objects
#'
#' Prediction helper for fitted GLBFP objects.
#'
#' @param object A fitted object of class `"glbfp_fit"` or `"glbfp_grid"`.
#' @param newdata Optional matrix/data frame with points where prediction is
#'   requested. For `"glbfp_fit"`, `newdata` is not supported.
#' @param ... Additional arguments (unused).
#'
#' @return Numeric vector of predicted densities.
#' @export
predict.glbfp_fit <- function(object, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    stop("`newdata` is not supported for single-point fit objects. Use *_estimate() for grid predictions.", call. = FALSE)
  }
  as.numeric(object$estimation)
}

#' @export
predict.glbfp_grid <- function(object, newdata = NULL, ...) {
  glbfp_nearest_grid_prediction(
    grid = object$grid,
    densities = object$densities,
    newdata = newdata
  )
}
