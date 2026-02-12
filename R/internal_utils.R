# Internal helpers for validation and grid handling.

#' @keywords internal
glbfp_validate_data <- function(data, allow_na = FALSE) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("`data` must be a matrix or data frame.", call. = FALSE)
  }

  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "double"

  if (nrow(data_mat) < 2L) {
    stop("`data` must contain at least two rows.", call. = FALSE)
  }
  if (ncol(data_mat) < 1L) {
    stop("`data` must contain at least one column.", call. = FALSE)
  }

  if (!allow_na && anyNA(data_mat)) {
    stop("`data` contains missing values. Remove or impute NA values first.", call. = FALSE)
  }
  if (any(!is.finite(data_mat))) {
    stop("`data` must contain only finite numeric values.", call. = FALSE)
  }

  data_mat
}

#' @keywords internal
glbfp_validate_vector <- function(x, d, name, positive = FALSE, integerish = FALSE) {
  if (!is.numeric(x) || length(x) != d || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a numeric vector of length %d with finite values.", name, d), call. = FALSE)
  }
  if (positive && any(x <= 0)) {
    stop(sprintf("`%s` must contain strictly positive values.", name), call. = FALSE)
  }
  if (integerish && any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop(sprintf("`%s` must contain integer values.", name), call. = FALSE)
  }

  if (integerish) {
    return(as.integer(round(x)))
  }
  as.numeric(x)
}

#' @keywords internal
glbfp_validate_bounds <- function(min_vals, max_vals, d) {
  min_vals <- glbfp_validate_vector(min_vals, d = d, name = "min_vals")
  max_vals <- glbfp_validate_vector(max_vals, d = d, name = "max_vals")
  if (any(min_vals >= max_vals)) {
    stop("Each entry in `min_vals` must be strictly smaller than `max_vals`.", call. = FALSE)
  }
  list(min_vals = min_vals, max_vals = max_vals)
}

#' @keywords internal
glbfp_cell_count <- function(lower, upper, step) {
  pmax(1L, as.integer(length(seq(lower, upper, by = step)) - 1L))
}

#' @keywords internal
glbfp_safe_estimation <- function(value) {
  value <- as.numeric(value)
  if (!is.finite(value) || is.na(value) || value < 0) {
    return(0)
  }
  value
}

#' @keywords internal
glbfp_prepare_grid <- function(data, grid_size = 20, grid_points = NULL, min_vals, max_vals) {
  d <- ncol(data)
  col_names <- colnames(data)
  if (is.null(col_names)) {
    col_names <- paste0("V", seq_len(d))
  }

  if (!is.null(grid_points)) {
    if (!is.matrix(grid_points) && !is.data.frame(grid_points)) {
      stop("`grid_points` must be a matrix or data frame when provided.", call. = FALSE)
    }
    grid <- as.matrix(grid_points)
    storage.mode(grid) <- "double"

    if (ncol(grid) != d) {
      stop("`grid_points` must have the same number of columns as `data`.", call. = FALSE)
    }
    if (any(!is.finite(grid))) {
      stop("`grid_points` must only contain finite values.", call. = FALSE)
    }

    grid_dims <- vapply(seq_len(d), function(i) length(unique(grid[, i])), integer(1))
    is_rectangular <- prod(grid_dims) == nrow(grid)
  } else {
    grid_size <- as.integer(round(grid_size))
    if (!is.finite(grid_size) || grid_size <= 1L) {
      stop("`grid_size` must be an integer greater than 1.", call. = FALSE)
    }

    axes <- lapply(seq_len(d), function(i) {
      seq(min_vals[i], max_vals[i], length.out = grid_size)
    })
    grid <- as.matrix(expand.grid(axes))
    grid_dims <- rep(grid_size, d)
    is_rectangular <- TRUE
  }

  colnames(grid) <- col_names

  list(
    grid = grid,
    grid_dims = grid_dims,
    is_rectangular = is_rectangular,
    col_names = col_names
  )
}

#' @keywords internal
glbfp_stabilize_covariance <- function(Sigma, max_attempts = 7L) {
  Sigma <- as.matrix(Sigma)
  d <- ncol(Sigma)
  Sigma <- 0.5 * (Sigma + t(Sigma))

  if (any(!is.finite(Sigma))) {
    stop("Covariance matrix contains non-finite values.", call. = FALSE)
  }

  diag_mean <- mean(diag(Sigma))
  if (!is.finite(diag_mean) || diag_mean <= 0) {
    diag_mean <- 1
  }

  ridge <- 0
  for (i in seq_len(max_attempts)) {
    Sigma_try <- Sigma + diag(ridge, d)
    Sigma_inv <- tryCatch(solve(Sigma_try), error = function(e) NULL)
    Sigma_det <- tryCatch(det(Sigma_try), error = function(e) NA_real_)

    if (!is.null(Sigma_inv) && is.finite(Sigma_det) && Sigma_det > 0 && all(is.finite(Sigma_inv))) {
      return(list(Sigma = Sigma_try, Sigma_inv = Sigma_inv, Sigma_det = Sigma_det, ridge = ridge))
    }

    ridge <- if (ridge == 0) sqrt(.Machine$double.eps) * diag_mean else ridge * 10
  }

  stop(
    "Failed to invert covariance matrix. `data` may be degenerate; try removing collinearity or reducing dimensionality.",
    call. = FALSE
  )
}

#' @keywords internal
glbfp_nearest_grid_prediction <- function(grid, densities, newdata) {
  if (is.null(newdata)) {
    return(densities)
  }

  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
    stop("`newdata` must be a matrix or data frame.", call. = FALSE)
  }
  newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"

  if (ncol(newdata) != ncol(grid)) {
    stop("`newdata` must have the same number of columns as the fitted grid.", call. = FALSE)
  }
  if (any(!is.finite(newdata))) {
    stop("`newdata` must contain only finite numeric values.", call. = FALSE)
  }

  idx <- vapply(seq_len(nrow(newdata)), function(i) {
    point <- newdata[i, ]
    distances <- rowSums((grid - matrix(point, nrow = nrow(grid), ncol = ncol(grid), byrow = TRUE))^2)
    which.min(distances)
  }, integer(1))

  as.numeric(densities[idx])
}
