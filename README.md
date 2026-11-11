<img src="man/figures/logo.svg" align="right" width="160" alt="GLBFP hex logo" />

# GLBFP

[![R-CMD-check](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml?query=branch%3Amain)
[![CRAN status](https://www.r-pkg.org/badges/version/GLBFP)](https://CRAN.R-project.org/package=GLBFP)
[![License](https://img.shields.io/cran/l/GLBFP)](https://CRAN.R-project.org/package=GLBFP)

GLBFP is an R package for nonparametric density estimation. It provides
Averaged Shifted Histogram (ASH), Linear Blend Frequency Polygon (LBFP), and
General Linear Blend Frequency Polygon (GLBFP) estimators.

The package estimates densities at individual points or on a grid and provides
functions for visualization and indexed random fields.

## Installation

Install [GLBFP 0.5.2 from CRAN](https://CRAN.R-project.org/package=GLBFP) for
ASH, LBFP and GLBFP density estimation.

```r
install.packages("GLBFP")
```

Install [GLBFP 0.6.1 from GitHub](https://github.com/AurelienNicosiaULaval/GLBFP/releases/tag/v0.6.1)
for indexed fields, including fields with more than one index dimension, and
`compute_Di()`.

```r
install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/GLBFP", ref = "v0.6.1")
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
| [`compute_Di()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_Di.html) | Compute fixed-grid leave-one-out diagnostics |

`glbfp_field()` and `compute_Di()` are available in the GitHub release.
For the assumptions behind the inference options, see `?LBFP`. The
`compute_Di()` score is a diagnostic, not an outlier probability. The LBFP
case is described in [Nicosia, Duchesne and Carbon (2026)](https://arxiv.org/abs/2607.14335).

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
- Nicosia, A., Duchesne, T. and Carbon, M. (2026). A Leave-One-Out Influence
  Statistic for Density-Based Outlier Detection. [arXiv:2607.14335](https://arxiv.org/abs/2607.14335).
- Scott, D. W. (1992). *Multivariate Density Estimation: Theory, Practice, and
  Visualization*. Wiley. [doi:10.1002/9780470316849](https://doi.org/10.1002/9780470316849).

Use `citation("GLBFP")` to cite the package.
