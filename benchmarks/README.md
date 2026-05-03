# Reproducible benchmarks

This directory contains benchmark scripts for development and manuscript
preparation. They are excluded from `R CMD check` through `.Rbuildignore`.

Run from the package root after installing development dependencies:

```r
install.packages(c("bench", "MASS"))
```

Then execute:

```sh
Rscript benchmarks/benchmark_glbfp_1d.R
Rscript benchmarks/benchmark_glbfp_2d.R
Rscript benchmarks/benchmark_loo_if_applicable.R
```

For a quick smoke test:

```sh
GLBFP_BENCH_QUICK=true Rscript benchmarks/benchmark_glbfp_1d.R
GLBFP_BENCH_QUICK=true Rscript benchmarks/benchmark_glbfp_2d.R
```

The scripts use fixed seeds and simulated data. They are intended to provide a
reproducible starting point for tables and figures in a software article. They
are not part of the package test suite.
