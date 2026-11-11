<img src="man/figures/logo.svg" align="right" width="160" alt="GLBFP hex logo" />

# GLBFP

[![R-CMD-check](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml/badge.svg?branch=develop)](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml?query=branch%3Adevelop)
[![CRAN status](https://www.r-pkg.org/badges/version/GLBFP)](https://CRAN.R-project.org/package=GLBFP)
[![License](https://img.shields.io/cran/l/GLBFP)](https://CRAN.R-project.org/package=GLBFP)

GLBFP is an R package for nonparametric density estimation. It provides
Averaged Shifted Histogram (ASH), Linear Blend Frequency Polygon (LBFP), and
General Linear Blend Frequency Polygon (GLBFP) estimators.

The package estimates densities at individual points or on a grid and provides
functions for visualization and indexed random fields.

## Installation

Install the stable version from CRAN:

```r
install.packages("GLBFP")
```

To use the development version, including indexed fields and `compute_Di()`:

```r
install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/GLBFP", ref = "develop")
```

## Example

```r
library(GLBFP)

set.seed(2026)
observations <- matrix(rnorm(300), ncol = 1)

fit <- GLBFP_estimate(data = observations, b = 0.5, m = 2, grid_size = 100)
plot(fit)
```

## Main functions

| Function | Purpose |
| --- | --- |
| `ASH()`, `LBFP()`, `GLBFP()` | Estimate a density at a point |
| `ASH_estimate()`, `LBFP_estimate()`, `GLBFP_estimate()` | Estimate a density on a grid |
| `compute_bi_optim()` | Select a reference bandwidth |
| `glbfp_field()` | Represent an indexed field |
| `plot_density()` | Display a density |
| `compute_Di()` | Compute fixed-grid leave-one-out diagnostics |

`glbfp_field()` and `compute_Di()` are available in the development version.
For interval estimation and its assumptions, see `?LBFP` in the installed
package.

## Documentation

Visit the [package website](https://aureliennicosiaulaval.github.io/GLBFP/)
for function documentation and examples, including
[Quebec climate fields](https://aureliennicosiaulaval.github.io/GLBFP/articles/quebec-climate-fields.html)
indexed in two and three dimensions.

## References

- Carbon, M. and Duchesne, T. (2024). Multivariate frequency polygon for
  stationary random fields. [doi:10.1007/s10463-023-00883-5](https://doi.org/10.1007/s10463-023-00883-5).
- Carbon, M. and Duchesne, T. (2026). Asymptotic normality of multivariate
  frequency polygons for stationary random fields. [doi:10.1007/s10463-025-00952-x](https://doi.org/10.1007/s10463-025-00952-x).

Use `citation("GLBFP")` to cite the package.
