# Fixed-seed experiment: independent, finite-range MA and Gaussian AR fields.
# Run from the package root after installing GLBFP 0.5.2.9002.
library(GLBFP)

set.seed(20261111)
replications <- 500L
bandwidth <- 0.3
anchor <- -5
query <- 0.15
truth <- dnorm(query)
shapes <- list(`1` = 900L, `2` = c(30L, 30L), `3` = c(10L, 10L, 10L))
rows <- vector("list", 3L * length(shapes) * replications)
row_id <- 0L
separable_ar <- function(shape, rho) {
  field <- array(rnorm(prod(shape)), dim = shape)
  for (axis in seq_along(shape)) {
    permutation <- c(axis, setdiff(seq_along(shape), axis))
    matrix_axis <- matrix(aperm(field, permutation), nrow = shape[axis])
    for (i in 2:shape[axis]) {
      matrix_axis[i, ] <- rho * matrix_axis[i - 1, ] + sqrt(1 - rho^2) * matrix_axis[i, ]
    }
    field <- aperm(array(matrix_axis, dim = shape[permutation]), order(permutation))
  }
  as.vector(field)
}
for (N in seq_along(shapes)) {
  shape <- shapes[[N]]
  sites <- as.matrix(expand.grid(lapply(shape, seq_len)))
  n <- nrow(sites)
  expanded <- shape + 1L
  strides <- c(1, head(cumprod(expanded), -1L))
  base <- 1L + as.vector((sites * 1) %*% strides)
  for (scenario in c("iid", "finite_ma", "strong_ar")) {
    rho <- switch(scenario, iid=0, finite_ma=.5, strong_ar=.9)
    for (r in seq_len(replications)) {
      if (scenario == "strong_ar") {
        values <- separable_ar(shape, rho)
      } else {
        innovations <- rnorm(prod(expanded))
        values <- innovations[base]
        if (rho != 0) {
          for (j in seq_len(N)) values <- values + rho * innovations[base - strides[j]]
        }
        values <- values / sqrt(1 + N * rho^2)
      }
      field <- glbfp_field(values, indices = sites,
                           dependence = if (rho == 0) "iid" else "dependent")
      fit <- withCallingHandlers(
        LBFP(query, field, b = bandwidth, min_vals = anchor, max_vals = 5,
             inference = "stationary_mixing"),
        warning = function(w) {
          if (startsWith(conditionMessage(w), "Stationary-mixing intervals are asymptotic")) invokeRestart("muffleWarning")
        }
      )
      stopifnot(all(is.finite(c(fit$estimation, fit$sd, fit$IC))))
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(N=N,n=n,scenario=scenario,rho=rho,replication=r,
        estimate=fit$estimation,standard_error=fit$sd,
        lower=fit$IC[1],upper=fit$IC[2],truth=truth,
        covered=fit$IC[1] <= truth && truth <= fit$IC[2])
    }
  }
}
raw <- do.call(rbind, rows)
summary <- do.call(rbind, lapply(split(raw, list(raw$N,raw$scenario), drop=TRUE),function(z) {
  coverage <- mean(z$covered)
  data.frame(N=z$N[1], n=z$n[1], scenario=z$scenario[1], rho=z$rho[1], replications=nrow(z),
    coverage=coverage, coverage_mc_se=sqrt(coverage*(1-coverage)/nrow(z)),
    mean_estimate=mean(z$estimate), truth=z$truth[1],
    empirical_sd=sd(z$estimate), mean_reported_se=mean(z$standard_error),
    mean_width=mean(z$upper-z$lower))
}))
rownames(summary) <- NULL
write.csv(raw,"benchmarks/results/2026-11-11/mixing-simulation-raw.csv",row.names=FALSE)
write.csv(summary,"benchmarks/results/2026-11-11/mixing-simulation-summary.csv",row.names=FALSE)
writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())),
           "benchmarks/results/2026-11-11/session-info.txt")
print(summary)
