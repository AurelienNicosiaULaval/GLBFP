# Compute the \\G^\*\\ bandwidth constant

Computes the dimension-dependent constant \\G^\*\\ used by
[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md).

## Usage

``` r
compute_G_star(d)
```

## Arguments

- d:

  A single positive integer giving the data dimension.

## Value

A positive numeric scalar.

## Details

The implemented formula is \$\$G^\* = 2^{\frac{3(d - 4)}{2(4 + d)}}
\exp\left(\frac{1}{4 + d}\right) \frac{\pi^{d / 2}}{4 + d}.\$\$

## See also

[`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md),
[`G_i()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/G_i.md),
[`K_mi()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/K_mi.md)

## Examples

``` r
compute_G_star(1)
#> [1] 0.2320261
compute_G_star(2)
#> [1] 0.4373872
```
