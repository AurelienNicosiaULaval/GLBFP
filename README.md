<img src="man/figures/logo.svg" align="right" width="160" alt="GLBFP hex logo" />

# GLBFP

[![R-CMD-check](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/AurelienNicosiaULaval/GLBFP/branch/main/graph/badge.svg)](https://app.codecov.io/gh/AurelienNicosiaULaval/GLBFP)
[![DOI](https://zenodo.org/badge/936140873.svg)](https://doi.org/10.5281/zenodo.17945962)

`GLBFP` is an R package for histogram-based nonparametric density
estimation. It implements:

- Averaged Shifted Histogram estimators: `ASH()`, `ASH_estimate()`
- Linear Blend Frequency Polygon estimators: `LBFP()`, `LBFP_estimate()`
- General Linear Blend Frequency Polygon estimators: `GLBFP()`, `GLBFP_estimate()`

The package supports pointwise density estimation, regular-grid estimation,
1D and 2D plotting, sparse-prefix grid-count computation, fixed-grid
leave-one-out `D_i` scores, S3 summaries and predictions, and plug-in bandwidth
selection.

## Status

The package is under active development and is not currently on CRAN. The
current development target is CRAN readiness and preparation of a software
article for The R Journal.

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

## Quick start

```r
library(GLBFP)

set.seed(2026)
x <- matrix(rnorm(300), ncol = 1)
b <- compute_bi_optim(x, m = 1)

fit <- glbfp(x = 0, data = x, b = b, m = 1)
fit
summary(fit)
predict(fit)
```

Uppercase function names remain available and are the historical API:

```r
fit_upper <- GLBFP(x = 0, data = x, b = b, m = 1)
identical(fit$estimation, fit_upper$estimation)
```

## Two-dimensional workflow

```r
library(GLBFP)

data("ashua")

river_data <- ashua[, c("flow", "level")]
b <- c(8, 0.4)
x0 <- c(mean(river_data$flow), mean(river_data$level))

point_fit <- glbfp(x = x0, data = river_data, b = b, m = c(1, 1))
point_fit

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

## Leave-one-out diagnostics

```r
scores <- compute_di(river_data, b = b, m = c(1, 1), estimator = "GLBFP")
summary(scores)
head(as.data.frame(scores))
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
| Bandwidth constants | `K_mi()`, `G_i()`, `compute_G_star()` |
| S3 helpers | `print()`, `summary()`, `predict()`, `plot()`, `as.data.frame()` |

## Documentation

The pkgdown site is organized as a reading path:

1. Getting started with GLBFP
2. Package overview and workflow map
3. Brief methodological background
4. Choosing between ASH, LBFP and GLBFP
5. Two-dimensional density estimation
6. Sparse-prefix computation
7. Leave-one-out `D_i` diagnostics
8. Objects, summaries and plotting
9. Validation and comparison
10. Legacy estimation example

The first five articles introduce the package and the estimators. The next
three articles document implementation diagnostics and S3 behavior. The
validation article gives a lightweight reproducible benchmark, while the legacy
vignette is kept for backward compatibility.

## References

General background on frequency polygons, averaged shifted histograms, and
multivariate density estimation is available in:

- Carbon, M., and Duchesne, T. (2024). Multivariate frequency polygon for
  stationary random fields. Annals of the Institute of Statistical Mathematics,
  76(2), 263-287. doi:10.1007/s10463-023-00883-5.
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
