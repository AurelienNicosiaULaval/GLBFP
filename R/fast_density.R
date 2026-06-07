# Fast sparse cell-count engines for ASH, LBFP and GLBFP.

#' @keywords internal
glbfp_key_strings <- function(keys) {
  keys <- as.matrix(keys)
  if (nrow(keys) == 0L) {
    return(character())
  }
  if (ncol(keys) == 1L) {
    return(as.character(keys[, 1]))
  }
  do.call(paste, c(as.data.frame(keys), sep = "\r"))
}

#' @keywords internal
glbfp_count_state <- function(data, delta, min_vals) {
  keys <- floor(sweep(sweep(data, 2, min_vals, "-"), 2, delta, "/")) + 1L
  keys <- as.matrix(keys)
  storage.mode(keys) <- "integer"

  key_values <- glbfp_key_strings(keys)
  counts <- table(key_values)
  count_keys <- names(counts)
  if (ncol(keys) == 1L) {
    cells <- matrix(as.integer(count_keys), ncol = 1L)
  } else {
    cells <- do.call(rbind, strsplit(count_keys, "\r", fixed = TRUE))
    storage.mode(cells) <- "integer"
  }
  count_values <- as.integer(counts)
  prefix_index <- glbfp_build_prefix_index(cells = cells, counts = count_values)

  list(
    counts = counts,
    cells = cells,
    count_values = count_values,
    prefix_index = prefix_index,
    data_keys = keys,
    delta = delta,
    min_vals = min_vals
  )
}

#' @keywords internal
glbfp_make_count_environment <- function(keys, counts) {
  env <- new.env(parent = emptyenv(), hash = TRUE, size = max(29L, length(keys) * 2L))
  for (i in seq_along(keys)) {
    assign(keys[[i]], counts[[i]], envir = env)
  }
  env
}

#' @keywords internal
glbfp_lookup_count_env <- function(env, key) {
  val <- mget(key, envir = env, ifnotfound = list(0L), inherits = FALSE)
  as.integer(unlist(val, use.names = FALSE))
}

#' @keywords internal
glbfp_choose_prefix_order <- function(cells) {
  n_unique <- apply(cells, 2L, function(z) length(unique(z)))
  order(n_unique, decreasing = FALSE)
}

#' @keywords internal
glbfp_build_prefix_index <- function(cells, counts, order = NULL) {
  cells <- as.matrix(cells)
  storage.mode(cells) <- "integer"
  d <- ncol(cells)
  if (is.null(order)) {
    order <- glbfp_choose_prefix_order(cells)
  }
  order <- as.integer(order)

  ordered_cells <- cells[, order, drop = FALSE]
  prefix_envs <- vector("list", d)
  for (depth in seq_len(d)) {
    prefix_keys <- glbfp_key_strings(ordered_cells[, seq_len(depth), drop = FALSE])
    env <- new.env(parent = emptyenv(), hash = TRUE, size = max(29L, length(prefix_keys) * 2L))
    for (key in unique(prefix_keys)) {
      assign(key, TRUE, envir = env)
    }
    prefix_envs[[depth]] <- env
  }

  full_keys <- glbfp_key_strings(ordered_cells)
  list(
    prefix_envs = prefix_envs,
    count_env = glbfp_make_count_environment(full_keys, counts),
    order = order
  )
}

#' @keywords internal
glbfp_prefix_exists <- function(prefix_index, key, depth) {
  exists(key, envir = prefix_index$prefix_envs[[depth]], inherits = FALSE)
}

#' @keywords internal
glbfp_sparse_prefix_eval <- function(state, cell_index, candidates, denominator, self_keys = NULL) {
  cell_index <- as.matrix(cell_index)
  storage.mode(cell_index) <- "integer"
  n_eval <- nrow(cell_index)
  d <- ncol(cell_index)
  prefix_index <- state$prefix_index
  order <- prefix_index$order

  densities <- numeric(n_eval)
  visited <- integer(n_eval)
  prefix_nodes <- integer(n_eval)
  self_weight <- if (is.null(self_keys)) NULL else numeric(n_eval)
  self_ordered_keys <- NULL
  if (!is.null(self_keys)) {
    self_ordered_keys <- glbfp_key_strings(self_keys[, order, drop = FALSE])
  }

  for (i in seq_len(n_eval)) {
    candidate_i <- if (is.function(candidates)) candidates(i) else candidates
    acc <- 0
    weight_self <- 0
    vi <- 0L
    nodes <- 0L

    rec <- function(depth, key, weight) {
      dim <- order[[depth]]
      cand <- candidate_i[[dim]]

      for (j in seq_along(cand$offset)) {
        w <- weight * cand$weight[[j]]
        if (w == 0) {
          next
        }

        coord <- cell_index[i, dim] + cand$offset[[j]]
        new_key <- if (depth == 1L) as.character(coord) else paste0(key, "\r", coord)
        nodes <<- nodes + 1L

        if (!glbfp_prefix_exists(prefix_index, new_key, depth)) {
          next
        }

        if (depth == d) {
          count <- glbfp_lookup_count_env(prefix_index$count_env, new_key)
          if (count != 0L) {
            vi <<- vi + 1L
            acc <<- acc + w * count
          }
          if (!is.null(self_ordered_keys) && identical(new_key, self_ordered_keys[[i]])) {
            weight_self <<- weight_self + w
          }
        } else {
          rec(depth + 1L, new_key, w)
        }
      }

      invisible(NULL)
    }

    rec(1L, "", 1)
    densities[[i]] <- acc / denominator
    visited[[i]] <- vi
    prefix_nodes[[i]] <- nodes
    if (!is.null(self_keys)) {
      self_weight[[i]] <- weight_self
    }
  }

  list(
    densities = densities,
    visited = visited,
    prefix_nodes = prefix_nodes,
    self_weight = self_weight,
    prefix_order = order
  )
}

#' @keywords internal
glbfp_evaluate_ash_fast <- function(data, points, b, m, min_vals, max_vals, self_keys = NULL) {
  points <- as.matrix(points)
  storage.mode(points) <- "double"

  d <- ncol(data)
  n <- nrow(data)
  delta <- b / m
  volume <- prod(b)
  state <- glbfp_count_state(data = data, delta = delta, min_vals = min_vals)

  cell_count <- vapply(seq_len(d), function(i) {
    glbfp_cell_count(min_vals[i], max_vals[i], delta[i])
  }, integer(1))

  n_eval <- nrow(points)
  cell_index <- matrix(NA_integer_, nrow = n_eval, ncol = d)
  for (i in seq_len(n_eval)) {
    cell_index[i, ] <- as.integer(pmin(
      pmax(1L, floor((points[i, ] - min_vals) / delta) + 1L),
      cell_count
    ))
  }

  candidates <- lapply(seq_len(d), function(s) {
    ell <- seq.int(1L - m[s], m[s] - 1L)
    list(offset = ell, weight = pmax(0, 1 - abs(ell) / m[s]))
  })
  sparse <- glbfp_sparse_prefix_eval(
    state = state,
    cell_index = cell_index,
    candidates = candidates,
    denominator = n * volume,
    self_keys = self_keys
  )

  list(
    densities = vapply(sparse$densities, glbfp_safe_estimation, numeric(1)),
    self_weight = sparse$self_weight,
    cell_index = cell_index,
    visited = sparse$visited,
    prefix_nodes = sparse$prefix_nodes,
    prefix_order = sparse$prefix_order,
    state = state
  )
}

#' @keywords internal
glbfp_evaluate_glbfp_fast <- function(data, points, b, m, min_vals, max_vals, self_keys = NULL) {
  points <- as.matrix(points)
  storage.mode(points) <- "double"

  d <- ncol(data)
  n <- nrow(data)
  delta <- b / m
  volume <- prod(b)
  state <- glbfp_count_state(data = data, delta = delta, min_vals = min_vals)

  a <- min_vals + delta / 2
  upper_vals <- max_vals + delta / 2
  cell_count <- vapply(seq_len(d), function(i) {
    glbfp_cell_count(a[i], upper_vals[i], delta[i])
  }, integer(1))

  n_eval <- nrow(points)
  u_mat <- matrix(NA_real_, nrow = n_eval, ncol = d)
  cell_index <- matrix(NA_integer_, nrow = n_eval, ncol = d)

  for (i in seq_len(n_eval)) {
    idx <- pmin(
      pmax(1L, floor((points[i, ] - a) / delta) + 1L),
      cell_count
    )
    idx <- as.integer(idx)
    cell_index[i, ] <- idx

    lowerbound_cell <- a + (idx - 1L) * delta
    u <- (points[i, ] - lowerbound_cell) / delta
    u_mat[i, ] <- u
  }

  candidates <- function(i) {
    lapply(seq_len(d), function(s) {
      ell <- seq.int(1L - m[s], m[s])
      tri_left <- pmax(0, 1 - abs(ell) / m[s])
      tri_right <- pmax(0, 1 - abs(ell - 1L) / m[s])
      list(
        offset = ell,
        weight = (1 - u_mat[i, s]) * tri_left + u_mat[i, s] * tri_right
      )
    })
  }
  sparse <- glbfp_sparse_prefix_eval(
    state = state,
    cell_index = cell_index,
    candidates = candidates,
    denominator = n * volume,
    self_keys = self_keys
  )

  densities <- vapply(sparse$densities, glbfp_safe_estimation, numeric(1))
  sigma_hat2 <- vapply(seq_len(n_eval), function(i) {
    glbfp_safe_estimation(
      (1 / (n * volume)) *
        prod((2 * m^2 + 1 - 6 * u_mat[i, ] * (1 - u_mat[i, ])) / (3 * m^2)) *
        densities[i]
    )
  }, numeric(1))

  list(
    densities = densities,
    sd = sqrt(sigma_hat2),
    IC = cbind(
      IC_lower = densities + stats::qnorm(0.025) * sqrt(sigma_hat2) / sqrt(n * volume),
      IC_upper = densities + stats::qnorm(0.975) * sqrt(sigma_hat2) / sqrt(n * volume)
    ),
    self_weight = sparse$self_weight,
    u = u_mat,
    cell_index = cell_index,
    visited = sparse$visited,
    prefix_nodes = sparse$prefix_nodes,
    prefix_order = sparse$prefix_order,
    state = state
  )
}

#' @keywords internal
glbfp_compute_di_from_density <- function(density, self_weight, n, volume) {
  D <- rep(NA_real_, length(density))
  ok <- is.finite(density) & density > 0 & is.finite(self_weight)
  D[ok] <- (self_weight[ok] / (volume * density[ok]) - 1) / (n - 1)
  D
}
