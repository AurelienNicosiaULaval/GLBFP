# GLBFP validation and comparison

This vignette provides a lightweight, fully reproducible benchmark for:

- 1D scenarios: S1 normal, S2 bimodal mixture, S3 heavy-tail (`t` with 4
  d.f.)
- 2D scenarios: S4 bimodal correlated mixture, S5 correlated Gaussian

Compared methods:

- [`ASH_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/ASH_estimate.md)
- [`LBFP_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/LBFP_estimate.md)
- [`GLBFP_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/GLBFP_estimate.md)
- KDE ([`stats::density`](https://rdrr.io/r/stats/density.html),
  [`MASS::kde2d`](https://rdrr.io/pkg/MASS/man/kde2d.html))

Metrics:

- ISE per run, MISE across repetitions
- elapsed time per run
- sensitivity to parameter changes (`b`, `m`)

## Run benchmark suite

``` r

result <- run_validation_suite(reps = 2, grid_size_1d = 60, grid_size_2d = 16)
summary_tbl <- summarise_validation(result)
summary_tbl
#>    scenario dimension method    ise_mean elapsed_mean  ise_median
#> 1        S1         1    ASH 0.032789580       0.0150 0.032789580
#> 2        S1         1  GLBFP 0.021953737       0.0245 0.021953737
#> 3        S1         1    KDE 0.002658040       0.0010 0.002658040
#> 4        S1         1   LBFP 0.021953737       0.0135 0.021953737
#> 5        S2         1    ASH 0.020182782       0.0140 0.020182782
#> 6        S2         1  GLBFP 0.013839910       0.0240 0.013839910
#> 7        S2         1    KDE 0.005633544       0.0005 0.005633544
#> 8        S2         1   LBFP 0.013808461       0.0135 0.013808461
#> 9        S3         1    ASH 0.029570127       0.0145 0.029570127
#> 10       S3         1  GLBFP 0.021946772       0.0245 0.021946772
#> 11       S3         1    KDE 0.006869287       0.0010 0.006869287
#> 12       S3         1   LBFP 0.021946772       0.0130 0.021946772
#> 13       S4         2    ASH 0.038106086       0.0840 0.038106086
#> 14       S4         2  GLBFP 0.018294805       0.1890 0.018294805
#> 15       S4         2    KDE 0.003425851       0.0010 0.003425851
#> 16       S4         2   LBFP 0.018285314       0.0940 0.018285314
#> 17       S5         2    ASH 0.047296014       0.0790 0.047296014
#> 18       S5         2  GLBFP 0.023920100       0.1915 0.023920100
#> 19       S5         2    KDE 0.004657292       0.0010 0.004657292
#> 20       S5         2   LBFP 0.023346380       0.0910 0.023346380
#>    elapsed_median
#> 1          0.0150
#> 2          0.0245
#> 3          0.0010
#> 4          0.0135
#> 5          0.0140
#> 6          0.0240
#> 7          0.0005
#> 8          0.0135
#> 9          0.0145
#> 10         0.0245
#> 11         0.0010
#> 12         0.0130
#> 13         0.0840
#> 14         0.1890
#> 15         0.0010
#> 16         0.0940
#> 17         0.0790
#> 18         0.1915
#> 19         0.0010
#> 20         0.0910
```

## MISE ranking by scenario

``` r

mi <- summary_tbl[order(summary_tbl$scenario, summary_tbl$ise_mean), ]
mi[, c("scenario", "dimension", "method", "ise_mean", "elapsed_mean")]
#>    scenario dimension method    ise_mean elapsed_mean
#> 3        S1         1    KDE 0.002658040       0.0010
#> 4        S1         1   LBFP 0.021953737       0.0135
#> 2        S1         1  GLBFP 0.021953737       0.0245
#> 1        S1         1    ASH 0.032789580       0.0150
#> 7        S2         1    KDE 0.005633544       0.0005
#> 8        S2         1   LBFP 0.013808461       0.0135
#> 6        S2         1  GLBFP 0.013839910       0.0240
#> 5        S2         1    ASH 0.020182782       0.0140
#> 11       S3         1    KDE 0.006869287       0.0010
#> 10       S3         1  GLBFP 0.021946772       0.0245
#> 12       S3         1   LBFP 0.021946772       0.0130
#> 9        S3         1    ASH 0.029570127       0.0145
#> 15       S4         2    KDE 0.003425851       0.0010
#> 16       S4         2   LBFP 0.018285314       0.0940
#> 14       S4         2  GLBFP 0.018294805       0.1890
#> 13       S4         2    ASH 0.038106086       0.0840
#> 19       S5         2    KDE 0.004657292       0.0010
#> 20       S5         2   LBFP 0.023346380       0.0910
#> 18       S5         2  GLBFP 0.023920100       0.1915
#> 17       S5         2    ASH 0.047296014       0.0790
```

## Sensitivity analysis (`b`, `m`)

``` r

result$sensitivity
#>   b_scale m mean_density max_density
#> 1     0.8 1   0.01875666   0.1965201
#> 2     1.0 1   0.01847158   0.1951697
#> 3     1.2 1   0.01785758   0.1301115
#> 4     1.0 2   0.01783902   0.1514042
```

## Notes

- The benchmark is intentionally small to keep vignette runtime short.
- Increase `reps` and grid sizes in `inst/bench/run_validation.R` for
  paper-grade tables.
