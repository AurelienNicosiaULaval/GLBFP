# GLBFP 0.5.1

## Features

- Added test infrastructure (`testthat` edition 3) and regression/invariant tests.
- Added S3 helpers for `summary()` and `predict()` on fit/grid objects.
- Added two publication-oriented vignettes (introduction and validation).
- Added reproducible benchmark scripts in `inst/bench`.
- Added pkgdown and coverage workflows.
- Added citation files (`CITATION.cff`, `inst/CITATION`) and `cran-comments.md`.

## Fixes

- Strengthened input validation (finite numeric checks, shape checks, explicit errors).
- Improved numerical robustness in `compute_bi_optim()` for near-singular covariance matrices.
- Fixed plotting methods for irregular custom grids and metadata handling.
- Fixed broken references in `plot.ASH_estimate()`.

# GLBFP 0.4.0

- Added standard deviation and confidence interval output for `*_estimate` functions.

# GLBFP 0.3.0

- Internal maintenance release.

# GLBFP 0.2.0

- Fixed errors in `ASH()` and `GLBFP()`.

# GLBFP 0.1.0

- Initial release.
