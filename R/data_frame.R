#' Convert GLBFP objects to data frames
#'
#' @param x A GLBFP grid object or a `compute_Di()` result.
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

#' @rdname as.data.frame.glbfp_grid
#' @export
as.data.frame.glbfp_di <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- data.frame(
    observation = seq_len(x$n),
    D = as.numeric(x$D),
    D_positive = as.numeric(x$D_positive),
    density = as.numeric(x$density),
    density_loo = as.numeric(x$density_loo),
    self_weight = as.numeric(x$self_weight),
    visited = as.integer(x$visited),
    prefix_nodes = as.integer(x$prefix_nodes)
  )
  row.names(out) <- row.names
  out
}
