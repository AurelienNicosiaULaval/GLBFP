#' Convert GLBFP objects to data frames
#'
#' @param x A GLBFP grid object.
#' @param row.names Optional row names.
#' @param optional Passed to [base::as.data.frame()].
#' @param ... Additional arguments (unused).
#'
#' @return A data frame representation of the object.
#'
#' @export
as.data.frame.glbfp_grid <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- as.data.frame(x$grid, row.names = row.names, optional = optional)
  names(out) <- x$col_names
  out$density <- as.numeric(x$densities)

  if (!is.null(x$sd)) {
    out$sd <- as.numeric(x$sd)
  }
  if (!is.null(x$IC)) {
    out$IC_lower <- as.numeric(x$IC[, "IC_lower"])
    out$IC_upper <- as.numeric(x$IC[, "IC_upper"])
  }
  if (!is.null(x$visited)) {
    out$visited <- as.integer(x$visited)
  }
  if (!is.null(x$prefix_nodes)) {
    out$prefix_nodes <- as.integer(x$prefix_nodes)
  }

  out
}
