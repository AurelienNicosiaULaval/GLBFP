# GLBFP estimation (legacy vignette)

This legacy vignette is kept for backward compatibility. For new users,
see Introduction to GLBFP and GLBFP validation and comparison.

``` r

library(GLBFP)
data("ashua")

x <- c(200, 30)
b <- c(0.5, 0.5)

fit <- GLBFP(x, ashua[, -3], b = b, m = c(1, 1))
fit
#> GLBFP Density Estimation:
#> Point: (200, 30) 
#> Estimated density: 0.003444976 
#> Estimated standard error: 0.00107138 
#> 95% confidence interval: 0.00338158345673361, 0.00350836869637647 
#> Bandwidths (b): 0.5, 0.5 
#> Shifts (m): 1, 1

grid_fit <- GLBFP_estimate(ashua[, -3], b = b, m = c(1, 1), grid_size = 12)
head(cbind(grid_fit$grid, density = grid_fit$densities))
#>       flow level density
#> [1,]  39.5 29.33       0
#> [2,] 175.0 29.33       0
#> [3,] 310.5 29.33       0
#> [4,] 446.0 29.33       0
#> [5,] 581.5 29.33       0
#> [6,] 717.0 29.33       0
```
