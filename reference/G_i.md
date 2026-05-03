# Compute the \\G(m_i)\\ bandwidth constant

Computes the scalar constant \\G(m_i)\\ used by
[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md).

## Usage

``` r
G_i(mi)
```

## Arguments

- mi:

  A single positive numeric value. In package estimators, `mi`
  corresponds to one component of the integer shift vector `m`.

## Value

A positive numeric scalar.

## Details

The implemented formula is \$\$G(m_i) = \frac{1}{12}\left(1 +
\frac{1}{2m_i^2}\right).\$\$

## See also

[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md),
[`K_mi()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/K_mi.md),
[`compute_G_star()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_G_star.md)

## Examples

``` r
G_i(1)
#> [1] 0.125
G_i(2)
#> [1] 0.09375
```
