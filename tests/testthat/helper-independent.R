# Independent observation-wise basis construction. No prefix traversal or
# evaluation-cell interpolation is used for the continuous estimators.
reference_weights <- function(data, points, b, m, origin, estimator = "GLBFP") {
  data <- as.matrix(data)
  points <- as.matrix(points)
  d <- ncol(data)
  if (estimator == "LBFP") m <- rep(1L, d)
  delta <- b / m
  keys <- floor(sweep(sweep(data, 2, origin, "-"), 2, delta, "/"))
  W <- matrix(1, nrow(points), nrow(data))
  for (j in seq_len(d)) {
    if (estimator == "ASH") {
      query_keys <- floor((points[, j] - origin[j]) / delta[j])
      w <- pmax(0, 1 - abs(outer(query_keys, keys[, j], "-")) / m[j])
    } else {
      centers <- origin[j] + (keys[, j] + 0.5) * delta[j]
      w <- matrix(0, nrow(points), nrow(data))
      for (shift in seq.int(1L - m[j], m[j] - 1L)) {
        w <- w + (1 - abs(shift) / m[j]) *
          pmax(0, 1 - abs(outer(points[, j], centers + shift * delta[j], "-")) / delta[j])
      }
    }
    W <- W * w
  }
  W / prod(b)
}

reference_density <- function(data, points, b, m, origin, estimator = "GLBFP") {
  rowMeans(reference_weights(data, points, b, m, origin, estimator))
}

reference_iid_variance <- function(W) apply(W, 1, stats::var) / ncol(W)

write_result <- function(x, file) utils::write.csv(x, file, row.names = FALSE)

timed <- function(fun, repetitions = 3L, minimum = 0.05) {
  # Warm up once, then repeat sufficiently to avoid sub-millisecond zero times.
  fun()
  vapply(seq_len(repetitions), function(i) {
    k <- 1L
    repeat {
      gc(FALSE)
      t <- system.time(for (j in seq_len(k)) fun())[["elapsed"]]
      if (t >= minimum || k >= 128L) return(t / k)
      k <- k * 2L
    }
  }, numeric(1))
}
