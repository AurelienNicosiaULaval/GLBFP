# Check the simulation generator against its known separable covariance.
set.seed(20261112)
separable_ar <- function(shape, rho) {
  field <- array(rnorm(prod(shape)), dim = shape)
  for (axis in seq_along(shape)) {
    permutation <- c(axis, setdiff(seq_along(shape), axis))
    matrix_axis <- matrix(aperm(field, permutation), nrow = shape[axis])
    for (i in 2:shape[axis])
      matrix_axis[i, ] <- rho * matrix_axis[i - 1, ] + sqrt(1 - rho^2) * matrix_axis[i, ]
    field <- aperm(array(matrix_axis, dim = shape[permutation]), order(permutation))
  }
  as.vector(field)
}
results <- lapply(1:3, function(N) {
  shape <- rep(3L,N)
  sites <- as.matrix(expand.grid(rep(list(1:3),N)))
  simulated <- replicate(5000L,separable_ar(shape,.9))
  empirical <- cov(t(simulated))
  expected <- outer(seq_len(nrow(sites)),seq_len(nrow(sites)),Vectorize(function(i,j) .9^sum(abs(sites[i,]-sites[j,]))))
  error <- max(abs(empirical-expected))
  stopifnot(error < .1)
  data.frame(N=N,replications=5000,max_absolute_error=error)
})
write.csv(do.call(rbind,results),"benchmarks/results/2026-11-11/generator-check.csv",row.names=FALSE)
print(do.call(rbind,results))
