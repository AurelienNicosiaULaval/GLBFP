# General Linear Blend Frequency Polygon (GLBFP) estimator at a single point

Computes the GLBFP density estimate at point `x`.

## Usage

``` r
GLBFP(
  x,
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
)

# S3 method for class 'GLBFP'
print(x, ...)
```

## Arguments

- x:

  Object returned by `GLBFP()`.

- data:

  Numeric matrix or data frame of observations (`n x d`).

- b:

  Positive numeric vector of bandwidths (length `d`).

- m:

  Positive integer vector of shifts (length `d`).

- min_vals:

  Numeric vector of lower grid bounds (length `d`).

- max_vals:

  Numeric vector of upper grid bounds (length `d`).

- ...:

  Additional arguments (unused).

## Value

A list with class `c("glbfp_fit", "GLBFP")` containing: `x`,
`estimation`, `sd`, `IC`, `b`, `m`, `method`, and `dimension`.

## Details

`GLBFP()` generalizes the linear blend frequency polygon workflow
through the positive integer shift vector `m`. Missing and non-finite
values are not accepted; remove or impute them before calling the
estimator.

## Methods (by generic)

- `print(GLBFP)`: Print method for object of class `"GLBFP"`.

## References

Scott, D. W. (1992). *Multivariate Density Estimation: Theory, Practice,
and Visualization*. Wiley. doi:10.1002/9780470316849.

The complete methodological citation for GLBFP has not yet been verified
in this repository. Add it before using this help page as publication
text.

## See also

[`GLBFP_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/GLBFP_estimate.md),
[`ASH()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/ASH.md),
[`LBFP()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/LBFP.md),
[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md)

## Examples

``` r
x <- c(200, 30)
b <- c(0.5, 0.5)
m <- c(1, 1)
GLBFP(x, ashua[, -3], b = b, m = m)
#> GLBFP Density Estimation:
#> Point: (200, 30) 
#> Estimated density: 0.003444976 
#> Estimated standard error: 0.00107138 
#> 95% confidence interval: 0.00338158345673361, 0.00350836869637647 
#> Bandwidths (b): 0.5, 0.5 
#> Shifts (m): 1, 1 
```
