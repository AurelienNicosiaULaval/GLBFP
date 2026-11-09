#' Prepare reusable sparse density counts
#'
#' Holds the fine-grid counts and an optional compiled index for repeated
#' evaluation. The R counts are retained so saved objects can rebuild their
#' native index after loading in another session.
#' @param data Numeric observations (`n x d`) or a [glbfp_field()].
#' @param b Positive bandwidth vector of length `d`. If `NULL`, use the
#'   Gaussian-reference L1 heuristic [compute_bi_optim()] with the supplied `m`.
#' @param m Positive integer vector of length `d`. `NULL` means all ones;
#'   the LBFP method always uses ones.
#' @param estimator `"GLBFP"`, `"LBFP"`, or `"ASH"`.
#' @param min_vals Fine-cell lower-bound anchor. `NULL` uses sample minima.
#' @param max_vals Upper plotting limits. They do not clamp the density support.
#'   `NULL` uses sample maxima. All cells outside the occupied set have zero count.
#' @param engine `"auto"`, `"R"`, `"cpp_direct"`, or `"cpp_prefix"`. Auto uses
#'   direct enumeration when the nominal stencil has at most 64 cells, prefix
#'   traversal otherwise. This is a computational heuristic, not an optimality
#'   claim. All engines use the same counts and weights.
#' @return A `glbfp_model` object, accepted by [compute_Di()] and [predict()].
#' @details
#' For `delta = b/m` and anchor `a = min_vals`, cell `k` (zero-based) is
#' `[a + k*delta, a + (k+1)*delta)`, with center `a + (k+1/2)*delta`.
#' ASH uses triangular weighted fine-cell counts. GLBFP interpolates ASH heights
#' at fine-cell centers. LBFP is GLBFP with `m=1`. The lattice is unbounded;
#' finite sample support is completed by zero counts, without extrapolation.
#' Exact boundary assignment uses floating-point `floor((x-a)/delta)` without
#' an undocumented snapping tolerance. Extremely small cells or unrepresentable
#' integer indices are rejected.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100), ncol = 1)
#' model <- prepare_glbfp(x, b = 0.8, m = 2)
#' predict(model, matrix(c(-1, 0, 1), ncol = 1))
#' compute_Di(model)
#' @export
prepare_glbfp <- function(data, b = NULL, m = NULL,
                          estimator = c("GLBFP", "LBFP", "ASH"),
                          min_vals = NULL, max_vals = NULL,
                          engine = c("auto", "R", "cpp_direct", "cpp_prefix")) {
  field <- if (inherits(data, "glbfp_field")) data else NULL
  data <- glbfp_validate_data(data)
  d <- ncol(data)
  estimator <- match.arg(estimator)
  engine <- match.arg(engine)
  fixed <- !is.null(b) && !is.null(min_vals)
  if (is.null(m) || estimator == "LBFP") m <- rep(1L, d)
  m <- glbfp_validate_vector(m, d, "m", positive = TRUE, integerish = TRUE)
  if (any(m > 1000000L)) stop("`m` is too large to enumerate a finite stencil safely.", call. = FALSE)
  if (is.null(b)) b <- compute_bi_optim(data, m = m)
  b <- glbfp_validate_vector(b, d, "b", positive = TRUE)
  delta <- b / m
  volume <- prod(b)
  if (any(!is.finite(delta) | delta <= 0) || !is.finite(volume) || volume <= 0 ||
      !is.finite(nrow(data) * volume)) {
    stop("Bandwidths produce a nonrepresentable cell width or volume.", call. = FALSE)
  }
  if (is.null(min_vals)) min_vals <- apply(data, 2, min)
  if (is.null(max_vals)) max_vals <- apply(data, 2, max)
  bounds <- glbfp_validate_bounds(min_vals, max_vals, d)
  # Positive delta alone does not ensure that cell centers can be represented
  # at the absolute scale of the data (for example 1e10 with delta=1e-12).
  locations <- rbind(bounds$min_vals, apply(data, 2, min), apply(data, 2, max))
  centers <- sweep(locations, 2, delta / 2, "+")
  next_edges <- sweep(locations, 2, delta, "+")
  if (any(!is.finite(centers)) || any(!is.finite(next_edges)) ||
      any(centers == locations) || any(next_edges == locations)) {
    stop("Cannot resolve cell centers at the observed scale; center/rescale data or increase bandwidths.", call. = FALSE)
  }
  nominal <- prod(if (estimator == "ASH") 2 * as.double(m) - 1 else 2 * as.double(m))
  if (engine == "auto") engine <- if (nominal <= 64) "cpp_direct" else "cpp_prefix"
  state <- glbfp_count_state(data, delta, bounds$min_vals, build_prefix = engine == "R")
  order <- glbfp_choose_prefix_order(state$cells)
  pointer <- if (engine == "R") NULL else glbfp_index_cpp(state$cells, state$count_values,
                                 as.integer(order), engine == "cpp_prefix")
  structure(list(data = data, state = state, pointer = pointer, engine = engine,
    b = b, m = m, delta = delta, volume = volume, min_vals = bounds$min_vals,
    max_vals = bounds$max_vals, n = nrow(data), d = d, method = estimator,
    prefix_order = order, nominal_stencil = nominal, field = field,
    fixed_parameters = fixed), class = "glbfp_model")
}

glbfp_evaluate_prepared <- function(model, points, self = FALSE, inference = "none") {
  points <- glbfp_validate_points(points, model$d)
  if (model$engine == "cpp_direct" && nrow(points) * model$nominal_stencil > 2e7) {
    stop("Direct stencil exceeds 20 million combinations; use `cpp_prefix` or smaller query batches.", call. = FALSE)
  }
  inference <- match.arg(inference, c("none", "iid"))
  if (inference == "iid") {
    if (!model$fixed_parameters) stop("IID uncertainty requires explicit fixed `b` and `min_vals`.", call. = FALSE)
    if (!is.null(model$field) && model$field$dependence != "iid") {
      stop("IID uncertainty requires `dependence = 'iid'` for an indexed field; dependent-field inference is not implemented.", call. = FALSE)
    }
  }
  ash <- model$method == "ASH"
  offset <- if (ash) 0 else model$delta / 2
  z <- sweep(sweep(points, 2, model$min_vals + offset, "-"), 2, model$delta, "/")
  outside <- rowSums(!is.finite(z) | abs(z) > .Machine$integer.max / 2) > 0
  z[outside, ] <- 0
  base <- floor(z) + 1L
  u <- z - floor(z)
  storage.mode(base) <- "integer"
  keys <- if (self) model$state$data_keys else matrix(integer(), nrow = 0, ncol = model$d)
  if (self && nrow(points) != model$n) stop("Self weights require the training observations.", call. = FALSE)
  if (model$engine == "R") {
    candidates <- function(i) lapply(seq_len(model$d), function(j) {
      ell <- seq.int(1L - model$m[j], if (ash) model$m[j] - 1L else model$m[j])
      left <- pmax(0, 1 - abs(ell) / model$m[j])
      list(offset = ell, weight = if (ash) left else
        (1 - u[i, j]) * left + u[i, j] * pmax(0, 1 - abs(ell - 1L) / model$m[j]))
    })
    out <- glbfp_sparse_prefix_eval(model$state, base, candidates, model$n * model$volume,
                                   if (self) keys else NULL)
  } else {
    if (!glbfp_index_valid_cpp(model$pointer)) {
      model$pointer <- glbfp_index_cpp(model$state$cells, model$state$count_values,
        as.integer(model$prefix_order), model$engine == "cpp_prefix")
    }
    out <- glbfp_eval_cpp(model$pointer, base, u, model$m, model$n, model$volume, ash, keys)
  }
  if (any(!is.finite(out$densities))) stop("Nonrepresentable density; rescale observations and bandwidths.", call. = FALSE)
  out$densities[outside] <- 0
  out$sum_square[outside] <- 0
  if (!is.null(out$self_weight)) out$self_weight[outside] <- 0
  # Unbiased sample estimate of Var(mean W) under IID sampling and fixed
  # parameters and evaluation points. The normal interval targets E[fhat],
  # with no claim to remove smoothing bias or guarantee finite-sample coverage.
  variance <- if (inference == "iid") pmax(0,
    (out$sum_square / model$volume / model$volume - model$n * out$densities^2) /
      (model$n * (model$n - 1))) else rep(NA_real_, nrow(points))
  if (inference == "iid" && any(!is.finite(variance))) stop("Nonrepresentable IID variance; rescale observations and bandwidths.", call. = FALSE)
  out$sd <- sqrt(variance)
  out$IC <- cbind(IC_lower = out$densities - stats::qnorm(.975) * out$sd,
                   IC_upper = out$densities + stats::qnorm(.975) * out$sd)
  out$u <- u
  out$cell_index <- base
  out$prefix_order <- model$prefix_order
  out$inference <- inference
  out
}

#' @rdname prepare_glbfp
#' @param object Prepared model.
#' @param newdata Matrix/data frame of evaluation points. `NULL` uses the
#'   training observations. Evaluation does not use nearest-grid interpolation.
#' @param ... Unused additional arguments.
#' @export
predict.glbfp_model <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data
  glbfp_evaluate_prepared(object, newdata)$densities
}

#' @export
print.glbfp_model <- function(x, ...) {
  cat(x$method, "prepared model:", x$n, "observations, d =", x$d, "\n")
  cat("Occupied fine cells:", nrow(x$state$cells), "; engine:", x$engine, "\n")
  invisible(x)
}
