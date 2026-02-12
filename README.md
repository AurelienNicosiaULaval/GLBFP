# GLBFP

[![R-CMD-check](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/AurelienNicosiaULaval/GLBFP/branch/main/graph/badge.svg)](https://app.codecov.io/gh/AurelienNicosiaULaval/GLBFP)
[![DOI](https://zenodo.org/badge/936140873.svg)](https://doi.org/10.5281/zenodo.17945962)

`GLBFP` is an R package for nonparametric density estimation using:

- ASH (`ASH`)
- LBFP (`LBFP`)
- GLBFP (`GLBFP`)

The package includes pointwise estimators, grid-based estimators, plotting
methods for 1D/2D, and a plug-in bandwidth selector.

## Installation

```r
install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/GLBFP")
```

## Quick start

```r
library(GLBFP)

data("ashua")
x <- c(200, 30)
b <- c(0.5, 0.5)

fit <- GLBFP(x = x, data = ashua[, -3], b = b, m = c(1, 1))
fit

grid_fit <- GLBFP_estimate(data = ashua[, -3], b = b, m = c(1, 1), grid_size = 20)
plot(grid_fit, contour = TRUE)
```

## Development checks

```r
devtools::document()
devtools::test()
devtools::check(manual = FALSE)
```

## Contributing

Please open an issue or pull request at:
<https://github.com/AurelienNicosiaULaval/GLBFP>
