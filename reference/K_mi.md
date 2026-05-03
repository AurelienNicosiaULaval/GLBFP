# Compute the \\K(m_i)\\ bandwidth constant

Computes the scalar constant \\K(m_i)\\ used by
[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md).

## Usage

``` r
K_mi(mi)
```

## Arguments

- mi:

  A single numeric value greater than 0.5. In package estimators, `mi`
  corresponds to one component of the integer shift vector `m`.

## Value

A positive numeric scalar.

## Details

The implemented formula is \$\$ K(m_i) = \sqrt{\frac{1}{6} +
\frac{1}{12m_i^2}} + \frac{4m_i^2 - 1}{6\sqrt{2}m_i} \log\left(
\frac{\sqrt{3} + \sqrt{4m_i^2 + 2}}{\sqrt{4m_i^2 - 1}} \right). \$\$

## See also

[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md),
[`G_i()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/G_i.md),
[`compute_G_star()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_G_star.md)

## Examples

``` r
K_mi(1)
#> [1] 0.8116126
K_mi(2)
#> [1] 0.8161827
```
