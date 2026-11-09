glbfp_validate_points <- function(points, d) {
  if (!is.matrix(points) && !is.data.frame(points)) stop("Evaluation points must be a matrix or data frame.", call. = FALSE)
  if (is.data.frame(points) && !all(vapply(points, is.numeric, logical(1)))) {
    stop("Evaluation points must contain numeric columns.", call. = FALSE)
  }
  points <- as.matrix(points)
  if (!is.numeric(points) || ncol(points) != d || any(!is.finite(points))) {
    stop("Evaluation points must be finite numeric values with d columns.", call. = FALSE)
  }
  storage.mode(points) <- "double"
  points
}

glbfp_fit_point <- function(x, data, b, m, method, min_vals, max_vals, engine, inference) {
  model <- prepare_glbfp(data, b, m, method, min_vals, max_vals, engine)
  x <- glbfp_validate_vector(x, model$d, "x")
  fast <- glbfp_evaluate_prepared(model, matrix(x, nrow = 1L), inference = inference)
  result <- list(x = x, estimation = fast$densities[1], sd = fast$sd[1], IC = fast$IC[1, ],
    b = model$b, m = model$m, method = method, dimension = model$d,
    u = if (method == "ASH") NULL else fast$u[1, ], cell_index = fast$cell_index[1, ],
    visited = fast$visited[1], prefix_nodes = fast$prefix_nodes[1],
    prefix_order = fast$prefix_order, engine = model$engine, inference = inference,
    min_vals = model$min_vals, max_vals = model$max_vals,
    index_dimension = if (is.null(model$field)) NULL else model$field$N)
  class(result) <- c("glbfp_fit", method)
  result
}

glbfp_fit_grid <- function(data, grid_size, grid_points, b, m, method,
                           min_vals, max_vals, engine, inference) {
  if (inference == "iid" && is.null(grid_points) && is.null(max_vals)) {
    stop("IID uncertainty on an automatic grid requires explicit fixed `max_vals`.", call. = FALSE)
  }
  model <- prepare_glbfp(data, b, m, method, min_vals, max_vals, engine)
  info <- glbfp_prepare_grid(model$data, grid_size, grid_points, model$min_vals, model$max_vals)
  fast <- glbfp_evaluate_prepared(model, info$grid, inference = inference)
  result <- list(grid = info$grid, densities = fast$densities, sd = fast$sd, IC = fast$IC,
    b = model$b, m = model$m, method = method, dimension = model$d,
    grid_size = if (length(unique(info$grid_dims)) == 1L) info$grid_dims[1] else NA_integer_,
    grid_dims = info$grid_dims, is_rectangular = info$is_rectangular, col_names = info$col_names,
    u = if (method == "ASH") NULL else fast$u, cell_index = fast$cell_index,
    visited = fast$visited, prefix_nodes = fast$prefix_nodes, prefix_order = fast$prefix_order,
    engine = model$engine, inference = inference, min_vals = model$min_vals,
    max_vals = model$max_vals, prepared = model,
    units = if (is.null(model$field)) NULL else model$field$units,
    index_dimension = if (is.null(model$field)) NULL else model$field$N)
  class(result) <- c("glbfp_grid", paste0(method, "_estimate"))
  result
}
