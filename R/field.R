#' Represent observations indexed by an integer lattice
#'
#' Stores the observation dimension `d` separately from the index dimension `N`.
#' Index coordinates are never appended to the observations used for density
#' estimation. The representation does not establish stationarity or mixing.
#'
#' @param values Numeric vector, matrix or data frame with explicit `indices`,
#'   or a numeric array without `indices`.
#' @param indices Numeric matrix or data frame (`n x N`) of distinct integer
#'   sites. Negative indices are allowed. With `NULL`, array axes define sites.
#' @param value_dim Optional axis of an array containing the observation
#'   components. All other axes are index dimensions. With `NULL`, array entries
#'   are scalar observations, including entries of a matrix.
#' @param na.action `"fail"` rejects missing observations; `"omit"` removes
#'   entire sites with at least one missing component and records their indices.
#'   Non-finite or missing indices and infinite observations are always rejected.
#' @param units Optional character vector of length `d` giving observation units.
#' @param index_units Optional character vector of length `N` describing index
#'   units. Physical coordinates on a noninteger lattice require a separate map.
#' @param dependence One of `"unspecified"`, `"dependent"`, or `"iid"`.
#'   `"iid"` is an explicit assumption by the caller, not a diagnostic.
#' @return A `glbfp_field` list with `values`, `indices`, `d`, `N`, bounding-box
#'   `shape`, `complete`, and `omitted` sites. Incomplete fields are accepted
#'   for descriptive marginal density estimation; no values are imputed.
#' @examples
#' set.seed(11)
#' field <- glbfp_field(array(rnorm(60), c(5, 4, 3)))
#' field
#' GLBFP(0, field, b = 1, m = 2, min_vals = -5, max_vals = 5)
#' sites <- as.matrix(expand.grid(i = 1:4, j = 1:5))
#' glbfp_field(cbind(a = rnorm(20), b = rnorm(20)), indices = sites)
#' @export
glbfp_field <- function(values, indices = NULL, value_dim = NULL,
                        na.action = c("fail", "omit"), units = NULL,
                        index_units = NULL,
                        dependence = c("unspecified", "dependent", "iid")) {
  na.action <- match.arg(na.action)
  dependence <- match.arg(dependence)
  if (is.null(indices)) {
    dims <- dim(values)
    if (!is.numeric(values) || is.null(dims) || any(dims < 1L)) {
      stop("Without `indices`, `values` must be a nonempty numeric array.", call. = FALSE)
    }
    axes <- seq_along(dims)
    if (!is.null(value_dim)) {
      if (length(value_dim) != 1L || !is.numeric(value_dim) ||
          !is.finite(value_dim) || value_dim != round(value_dim) ||
          !(value_dim %in% axes) || length(dims) < 2L) {
        stop("`value_dim` must identify one array axis, with at least one index axis remaining.", call. = FALSE)
      }
      axes <- axes[axes != value_dim]
      values <- matrix(aperm(values, c(axes, value_dim)), ncol = dims[value_dim])
    } else values <- matrix(values, ncol = 1L)
    indices <- as.matrix(expand.grid(lapply(dims[axes], seq_len)))
    colnames(indices) <- paste0("i", seq_along(axes))
  } else {
    if (!is.null(value_dim)) stop("`value_dim` applies only to array input without `indices`.", call. = FALSE)
    if (is.numeric(values) && is.null(dim(values))) values <- matrix(values, ncol = 1L)
  }
  if (!is.matrix(indices) && !is.data.frame(indices)) {
    stop("`indices` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (is.data.frame(indices) && !all(vapply(indices, is.numeric, logical(1)))) {
    stop("`indices` must contain only numeric columns.", call. = FALSE)
  }
  indices <- as.matrix(indices)
  if (!is.numeric(indices) || ncol(indices) < 1L || any(!is.finite(indices)) ||
      any(indices != round(indices)) || any(abs(indices) > 2^53 - 1)) {
    stop("`indices` must contain finite exactly represented integer sites.", call. = FALSE)
  }
  if (anyDuplicated(as.data.frame(indices))) stop("Duplicate index sites are not allowed.", call. = FALSE)
  if ((!is.matrix(values) && !is.data.frame(values)) ||
      (is.data.frame(values) && !all(vapply(values, is.numeric, logical(1))))) {
    stop("`values` must contain only numeric observation columns.", call. = FALSE)
  }
  values <- as.matrix(values)
  if (!is.numeric(values) || ncol(values) < 1L || nrow(values) != nrow(indices)) {
    stop("Observation and index row counts must agree, with at least one observation component.", call. = FALSE)
  }
  if (any(is.infinite(values))) stop("Infinite observations are not allowed.", call. = FALSE)
  missing <- !stats::complete.cases(values)
  if (any(missing) && na.action == "fail") stop("`values` contains missing observations.", call. = FALSE)
  omitted <- indices[missing, , drop = FALSE]
  values <- values[!missing, , drop = FALSE]
  indices <- indices[!missing, , drop = FALSE]
  if (nrow(values) < 2L) stop("At least two observed sites must remain.", call. = FALSE)
  d <- ncol(values); N <- ncol(indices)
  for (entry in c("units", "index_units")) {
    u <- get(entry)
    if (!is.null(u) && (!is.character(u) || anyNA(u) ||
                       length(u) != if (entry == "units") d else N)) {
      stop(paste0("`", entry, "` has invalid type or length."), call. = FALSE)
    }
  }
  shape <- apply(indices, 2, function(z) max(z) - min(z) + 1)
  storage.mode(values) <- "double"
  structure(list(values = values, indices = indices, d = d, N = N,
    shape = shape, complete = prod(shape) == nrow(indices), omitted = omitted,
    units = units, index_units = index_units, dependence = dependence), class = "glbfp_field")
}

#' @export
print.glbfp_field <- function(x, ...) {
  cat("Indexed field:", nrow(x$values), "observed sites; N =", x$N, "; d =", x$d, "\n")
  cat("Complete integer bounding box:", x$complete, "\n")
  cat("Omitted sites:", nrow(x$omitted), "; dependence:", x$dependence, "\n")
  invisible(x)
}
