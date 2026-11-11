# Muffle only the documented asymptotic warning; retain all unexpected warnings.
with_mixing_warning <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (startsWith(conditionMessage(w), "Stationary-mixing intervals are asymptotic")) invokeRestart("muffleWarning")
})

test_that("LBFP field intervals match independent corner-count calculations", {
  for (N in 1:3) {
    indices <- as.matrix(expand.grid(rep(list(0:2), N)))
    observations <- rep(c(-0.4, 0.2, 0.8, 1.2, 1.7, 2.1), length.out = nrow(indices))
    field <- glbfp_field(observations, indices = indices, dependence = "dependent")
    x <- 0.7; b <- 1; anchor <- -1
    # The two histogram centres around x are 0.5 and 1.5.
    counts <- c(sum(observations >= 0 & observations < 1),
                sum(observations >= 1 & observations < 2))
    weights <- c(0.8, 0.2)
    density <- sum(weights * counts) / length(observations)
    variance <- density * sum(weights^2) / length(observations)
    for (engine in c("R", "cpp_direct", "cpp_prefix")) {
      expect_warning(fit <- LBFP(x, field, b = b, min_vals = anchor, max_vals = 3,
                  engine = engine, inference = "stationary_mixing"), "asymptotic")
      expect_equal(fit$estimation, density, tolerance = 1e-12)
      expect_equal(fit$sd, sqrt(variance), tolerance = 1e-12)
      expect_equal(unname(fit$IC), density + c(-1, 1) * qnorm(.975) * sqrt(variance), tolerance = 1e-12)
      expect_match(paste(capture.output(print(fit)), collapse = " "), "interval for density")
      expect_match(paste(capture.output(print(summary(fit))), collapse = " "), "interval for density")
    }
  }
})

test_that("multivariate phase factors use the product of corner weights", {
  observations <- rbind(c(.1,.1), c(.2,1.1), c(1.1,.1), c(1.2,1.1), c(.3,.4), c(1.3,1.4))
  indices <- as.matrix(expand.grid(i=1:2,j=1:3))
  field <- glbfp_field(observations, indices, dependence="dependent")
  # x = (0.7,1.1), centres (0.5,1.5): u=(0.2,0.6).
  x <- c(.7,1.1); weights <- c(.8*.4,.2*.4,.8*.6,.2*.6)
  counts <- c(2,1,1,2)
  density <- sum(weights * counts) / nrow(observations)
  expected_variance <- density * sum(weights^2) / nrow(observations)
  fit <- with_mixing_warning(LBFP(x,field,b=c(1,1),min_vals=c(-1,-1),max_vals=c(3,3),inference="stationary_mixing"))
  expect_equal(fit$estimation,density,tolerance=1e-12)
  expect_equal(fit$sd,sqrt(expected_variance),tolerance=1e-12)
  fit2 <- with_mixing_warning(GLBFP(x,field,b=c(1,1),m=c(1,1),min_vals=c(-1,-1),max_vals=c(3,3),inference="stationary_mixing"))
  expect_equal(fit2$IC,fit$IC)
  grid <- with_mixing_warning(LBFP_estimate(field,b=c(1,1),min_vals=c(-1,-1),max_vals=c(3,3),
                         grid_points=matrix(x,nrow=1),inference="stationary_mixing"))
  expect_equal(unname(grid$IC[1,]),unname(fit$IC))
})

test_that("mixing mode rejects geometries and estimators outside its scope", {
  field <- glbfp_field(array(seq_len(9),c(3,3)),dependence="dependent")
  expect_error(LBFP(1, field, b=1, inference="stationary_mixing"), "explicit fixed")
  expect_error(LBFP(1, matrix(1:9,ncol=1), b=1, min_vals=0, inference="stationary_mixing"), "glbfp_field")
  unspecified <- glbfp_field(array(seq_len(9),c(3,3)))
  expect_error(LBFP(1,unspecified,b=1,min_vals=0,inference="stationary_mixing"), "Declare")
  incomplete <- glbfp_field(1:3,indices=matrix(c(1,1,2,2,3,1),ncol=2,byrow=TRUE),dependence="dependent")
  expect_error(LBFP(1,incomplete,b=1,min_vals=0,inference="stationary_mixing"), "complete rectangular")
  missing <- array(1:9,c(3,3)); missing[1,1] <- NA
  missing <- glbfp_field(missing,na.action="omit",dependence="dependent")
  expect_error(LBFP(1,missing,b=1,min_vals=0,inference="stationary_mixing"), "complete rectangular")
  expect_error(GLBFP(1,field,b=1,m=2,min_vals=0,inference="stationary_mixing"), "only for LBFP")
  expect_error(LBFP_estimate(field,b=1,min_vals=0,inference="stationary_mixing"), "fixed `max_vals`")
  expect_warning(fit <- with_mixing_warning(LBFP(100,field,b=1,min_vals=0,max_vals=10,inference="stationary_mixing")), "unavailable")
  expect_true(is.na(fit$sd))
  expect_true(all(is.na(fit$IC)))
})
