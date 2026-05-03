# Averaged Shifted Histogram (ASH) estimator at a single point

Computes the ASH density estimate at point `x`.

## Usage

``` r
ASH(
  x,
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
)

# S3 method for class 'ASH'
print(x, ...)
```

## Arguments

- x:

  Object of class `"ASH"`.

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

A list with class `c("glbfp_fit", "ASH")` containing: `x`, `estimation`,
`b`, `m`, `method`, and `dimension`.

## Details

`m` controls the number of shifted histograms used in each dimension.
Missing and non-finite values are not accepted; remove or impute them
before calling the estimator.

## Methods (by generic)

- `print(ASH)`: Print method for object of class `"ASH"`.

## References

Scott, D. W. (1992). *Multivariate Density Estimation: Theory, Practice,
and Visualization*. Wiley. doi:10.1002/9780470316849.

## See also

[`ASH_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/ASH_estimate.md),
[`LBFP()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/LBFP.md),
[`GLBFP()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/GLBFP.md),
[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md)

## Examples

``` r
x <- c(200, 30)
b <- c(0.5, 0.5)
m <- c(1, 1)
ASH(x, ashua[, -3], b = b, m = m)
#> ASH Density Estimation:
#> Point: (200, 30) 
#> Estimated density: 0.008202324 
#> Bandwidths (b): 0.5, 0.5 
#> Shifts (m): 1, 1 
```
