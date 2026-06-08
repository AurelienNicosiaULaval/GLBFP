## Test environments

- local macOS (aarch64, R 4.5.0)
- GitHub Actions: macOS-latest (release), windows-latest (release), ubuntu-latest (devel/release/oldrel-1)

## R CMD check results

`R CMD check --as-cran` was run on local source tarball.

- 0 ERROR
- 0 WARNING
- 1 NOTE

The NOTE is the expected CRAN incoming feasibility note for a new submission.

## Downstream and additional checks

- Test suite run file by file with `testthat::test_file()`
- Source archive built with `R CMD build` from a clean copy of the package
- `R CMD check --as-cran` run on the source archive
- URL check via `urlchecker::url_check()`

## Resubmission

Initial submission.
