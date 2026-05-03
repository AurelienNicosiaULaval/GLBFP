# Changelog

## GLBFP 0.5.1

### Features

- Added test infrastructure (`testthat` edition 3) and
  regression/invariant tests.
- Added S3 helpers for
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) on fit/grid
  objects.
- Added two publication-oriented vignettes (introduction and
  validation).
- Added reproducible benchmark scripts in `benchmarks/` and validation
  helpers in `inst/bench`.
- Added pkgdown and coverage workflows.
- Added citation files (`CITATION.cff`, `inst/CITATION`) and
  `cran-comments.md`.
- Added pre-R Journal preparation files in `dev/` and `rjournal-prep/`.

### Fixes

- Strengthened input validation (finite numeric checks, shape checks,
  explicit errors).
- Improved numerical robustness in
  [`compute_bi_optim()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/compute_bi_optim.md)
  for near-singular covariance matrices.
- Fixed plotting methods for irregular custom grids and metadata
  handling.
- Fixed broken references in
  [`plot.ASH_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/ASH_estimate.md).

## GLBFP 0.4.0

- Added standard deviation and confidence interval output for
  `*_estimate` functions.

## GLBFP 0.3.0

- Internal maintenance release.

## GLBFP 0.2.0

- Fixed errors in
  [`ASH()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/ASH.md)
  and
  [`GLBFP()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/GLBFP.md).

## GLBFP 0.1.0

- Initial release.
