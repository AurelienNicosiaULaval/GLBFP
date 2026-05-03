# Compute bandwidth vector \\b_i\\

Computes a plug-in bandwidth vector used by GLBFP/LBFP/ASH estimators.
The function validates numeric inputs, stabilizes near-singular
covariance matrices with a small ridge if needed, and returns strictly
positive bandwidths.

## Usage

``` r
compute_bi_optim(data, m = rep(1, ncol(data)))
```

## Arguments

- data:

  A numeric matrix or data frame where rows are observations and columns
  are variables.

- m:

  A positive integer vector of shifts, one value per dimension.

## Value

A numeric vector of positive bandwidths with one value per column in
`data`.

## Details

The returned vector is intended as a starting value for examples and
routine workflows. For applied analysis, sensitivity to the bandwidth
should still be checked.

Near-singular covariance matrices are stabilized with a small ridge
term. If this fails, the function returns an error rather than silently
producing non-finite bandwidths.

## See also

[`ASH()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/ASH.md),
[`LBFP()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/LBFP.md),
[`GLBFP()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/GLBFP.md)

## Examples

``` r
set.seed(1)
x <- cbind(rnorm(200), rnorm(200))
compute_bi_optim(x, m = c(1, 1))
#> [1] 0.3109093 0.3153027
```
