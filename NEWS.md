# GLBFP 0.5.2.9000

- Added `compute_Di()` for fixed-grid leave-one-out
  self-support scores. This function is available only in the GitHub
  development version and are not part of CRAN 0.5.2.
- Added data-frame, summary and plotting methods for `D_i` diagnostics, tests
  against explicit leave-one-out recomputation, and a dedicated vignette.
- Linked the single-grid LBFP diagnostic to the accompanying article by
  Nicosia, Duchesne and Carbon (2026), arXiv:2607.14335.
- Updated package authorship and citation metadata for the GitHub development
  version.

# GLBFP 0.5.1

## Features

- Added sparse-prefix computation for ASH, LBFP and GLBFP grid-count estimators.
- Added lowercase aliases (`ash()`, `lbfp()`, `glbfp()` and corresponding grid functions).
- Added `as.data.frame()` helpers for grid objects.
- Added a hexagonal package logo and expanded pkgdown article structure.
- Added test infrastructure (`testthat` edition 3) and regression/invariant tests.
- Added S3 helpers for `summary()` and `predict()` on fit/grid objects.
- Added two publication-oriented vignettes (introduction and validation).
- Added reproducible benchmark scripts in `benchmarks/` and validation helpers in `inst/bench`.
- Added pkgdown and coverage workflows.
- Added citation files (`CITATION.cff`, `inst/CITATION`) and `cran-comments.md`.
- Added pre-R Journal preparation files in `dev/` and `rjournal-prep/`.

## Fixes

- Strengthened input validation (finite numeric checks, shape checks, explicit errors).
- Improved numerical robustness in `compute_bi_optim()` for near-singular covariance matrices.
- Fixed plotting methods for irregular custom grids and metadata handling.
- Fixed broken references in `plot.ASH_estimate()`.
- Reduced duplicated S3 plotting and printing code through shared internal helpers.

# GLBFP 0.4.0

- Added standard deviation and confidence interval output for `*_estimate` functions.

# GLBFP 0.3.0

- Internal maintenance release.

# GLBFP 0.2.0

- Fixed errors in `ASH()` and `GLBFP()`.

# GLBFP 0.1.0

- Initial release.
