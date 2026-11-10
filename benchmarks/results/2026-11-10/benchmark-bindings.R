library(utils)
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(normalizePath(args[1]), .libPaths()))
library(GLBFP)
set.seed(20261110)
cases <- data.frame(n = c(10000L, 2000L, 300L), d = c(1L, 2L, 3L), m = c(2L, 2L, 4L))
measure <- function(fun) vapply(seq_len(3L), function(repetition) {
  k <- 1L
  repeat {
    elapsed <- system.time(for (j in seq_len(k)) fun())[['elapsed']]
    if (elapsed >= .05 || k >= 512L) return(elapsed / k)
    k <- k * 2L
  }
}, numeric(1))
rows <- list()
for (case in seq_len(nrow(cases))) {
  x <- matrix(runif(cases$n[case] * cases$d[case], -1, 1), ncol = cases$d[case])
  b <- rep(.8, ncol(x)); m <- rep(cases$m[case], ncol(x))
  anchor <- rep(-5, ncol(x)); upper <- rep(5, ncol(x))
  model <- prepare_glbfp(x, b, m, min_vals = anchor, max_vals = upper)
  q <- x[seq_len(min(nrow(x), 500)), , drop=FALSE]
  functions <- list(prepare = function() prepare_glbfp(x, b, m, min_vals=anchor, max_vals=upper),
    evaluate = function() predict(model,q),
    all_scores_prepared = function() compute_Di(model),
    all_scores_total = function() compute_Di(x,b,m,min_vals=anchor,max_vals=upper))
  for (phase in names(functions)) {
    values <- measure(functions[[phase]])
    rows[[length(rows)+1L]] <- data.frame(version=args[2], case=case, n=nrow(x), d=ncol(x), m=cases$m[case], phase=phase,
      repetition=seq_along(values), seconds=values)
  }
}
write.csv(do.call(rbind,rows),args[3],row.names=FALSE)
