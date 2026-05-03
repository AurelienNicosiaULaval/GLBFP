# Linear Blend Frequency Polygon (LBFP) estimator at a single point

Computes the LBFP density estimate at point `x`.

## Usage

``` r
LBFP(
  x,
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
)

# S3 method for class 'LBFP'
print(x, ...)
```

## Arguments

- x:

  Object from `LBFP()`.

- data:

  Numeric matrix or data frame of observations (`n x d`).

- b:

  Positive numeric vector of bandwidths (length `d`).

- min_vals:

  Numeric vector of lower grid bounds (length `d`).

- max_vals:

  Numeric vector of upper grid bounds (length `d`).

- ...:

  Additional arguments (unused).

## Value

A list with class `c("glbfp_fit", "LBFP")` containing: `x`,
`estimation`, `sd`, `IC`, `b`, `method`, and `dimension`.

## Details

The estimate is obtained by linear blending of neighboring histogram bin
heights. Missing and non-finite values are not accepted; remove or
impute them before calling the estimator.

## Methods (by generic)

- `print(LBFP)`: Print method for object of class `"LBFP"`.

## References

Scott, D. W. (1992). *Multivariate Density Estimation: Theory, Practice,
and Visualization*. Wiley. doi:10.1002/9780470316849.

Terrell, G. R., and Scott, D. W. (1985). Oversmoothed Nonparametric
Density Estimates. *Journal of the American Statistical Association*,
80(389), 209-214. doi:10.1080/01621459.1985.10477163.

## See also

[`LBFP_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/LBFP_estimate.md),
[`ASH()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/ASH.md),
[`GLBFP()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/GLBFP.md),
[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md)

## Examples

``` r
x <- c(200, 30)
b <- c(0.5, 0.5)
LBFP(x, ashua[, -3], b = b)
#> LBFP Density Estimation:
#> Point: (200, 30) 
#> Estimated density: 0.003444976 
#> Estimated standard error: 0.00107138 
#> 95% confidence interval: 0.00338158345673361, 0.00350836869637647 
#> Bandwidths (b): 0.5, 0.5 
```
