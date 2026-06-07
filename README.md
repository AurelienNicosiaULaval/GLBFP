<img src="man/figures/logo.svg" align="right" width="160" alt="GLBFP hex logo" />

# GLBFP

[![R-CMD-check](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/AurelienNicosiaULaval/GLBFP/branch/main/graph/badge.svg)](https://app.codecov.io/gh/AurelienNicosiaULaval/GLBFP)
[![DOI](https://zenodo.org/badge/936140873.svg)](https://doi.org/10.5281/zenodo.17945962)

`GLBFP` provides R implementations of histogram-based nonparametric density
estimators:

- Averaged Shifted Histogram: `ASH()`, `ASH_estimate()`
- Linear Blend Frequency Polygon: `LBFP()`, `LBFP_estimate()`
- General Linear Blend Frequency Polygon: `GLBFP()`, `GLBFP_estimate()`

The package includes pointwise estimators, grid estimators, sparse-prefix
grid-count computation, fixed-grid leave-one-out \(D_i\) scores, plotting
methods for 1D/2D grids, S3 helpers, and a plug-in bandwidth selector.

## Status

The package is under active development and is not currently on CRAN. The current
development target is CRAN readiness and preparation of a software article for
The R Journal.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/GLBFP")
```

If the package is later accepted on CRAN, installation will use:

```r
install.packages("GLBFP")
```

## Minimal example

```r
library(GLBFP)

set.seed(2026)
x <- matrix(rnorm(300), ncol = 1)
b <- compute_bi_optim(x, m = 1)

fit <- GLBFP(x = 0, data = x, b = b, m = 1)
fit
summary(fit)
predict(fit)
```

## Two-dimensional example

```r
library(GLBFP)

data("ashua")

river_data <- ashua[, c("flow", "level")]
b <- c(8, 0.4)
x0 <- c(mean(river_data$flow), mean(river_data$level))

fit <- glbfp(x = x0, data = river_data, b = b, m = c(1, 1))
fit

grid_fit <- glbfp_estimate(
  data = river_data,
  b = b,
  m = c(1, 1),
  grid_size = 20
)

summary(grid_fit)
head(as.data.frame(grid_fit))
plot(grid_fit, contour = TRUE)
```

## Leave-one-out scores

```r
scores <- compute_di(river_data, b = b, m = c(1, 1), estimator = "GLBFP")
summary(scores)
plot(scores)
```

## Main functions

| Task | Functions |
|---|---|
| Pointwise density estimation | `ASH()`, `LBFP()`, `GLBFP()` |
| Grid-based density estimation | `ASH_estimate()`, `LBFP_estimate()`, `GLBFP_estimate()` |
| Lowercase aliases | `ash()`, `lbfp()`, `glbfp()`, `ash_estimate()`, `lbfp_estimate()`, `glbfp_estimate()` |
| Leave-one-out diagnostics | `compute_Di()`, `compute_di()` |
| Bandwidth helper | `compute_bi_optim()` |
| Constants used by the bandwidth helper | `K_mi()`, `G_i()`, `compute_G_star()` |
| Basic S3 helpers | `print()`, `plot()`, `summary()`, `predict()`, `as.data.frame()` |

## Vignettes and pkgdown site

The pkgdown site is configured with reproducible guides for:

- getting started with GLBFP;
- a short methodological overview;
- two-dimensional density estimation;
- estimator choice;
- sparse-prefix computation;
- leave-one-out \(D_i\) scores;
- S3 objects and helper methods;
- validation and benchmark scenarios.

## References

General background on frequency polygons, averaged shifted histograms, and
multivariate density estimation is available in:

- Scott, D. W. (1992). Multivariate Density Estimation: Theory, Practice, and
  Visualization. Wiley. doi:10.1002/9780470316849.
- Terrell, G. R., and Scott, D. W. (1985). Oversmoothed Nonparametric Density
  Estimates. Journal of the American Statistical Association, 80(389), 209-214.
  doi:10.1080/01621459.1985.10477163.

The complete bibliographic record for the original GLBFP methodological article
has not yet been verified in this repository. It is tracked in
`dev/references_to_verify.md` and should be added before journal submission.

## Citation

To cite the package from R:

```r
citation("GLBFP")
```

The repository also includes `CITATION.cff` for software citation metadata.

## Development checks

```r
devtools::document()
devtools::test()
devtools::check()
rcmdcheck::rcmdcheck(args = c("--as-cran"))
```

Benchmarks are stored in `benchmarks/` and are not run automatically during
`R CMD check`.

## Contributing

Please use GitHub issues and pull requests:
<https://github.com/AurelienNicosiaULaval/GLBFP/issues>

## License

GPL (>= 3).
