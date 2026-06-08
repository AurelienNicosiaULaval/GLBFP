## Test environments

- local macOS (aarch64, R 4.5.0)
- GitHub Actions: macOS-latest (release), windows-latest (release), ubuntu-latest (devel/release/oldrel-1)

## R CMD check results

`R CMD check --as-cran` was run on local source tarball.

- 0 ERROR
- 0 WARNING
- 2 NOTEs

The first NOTE is the expected CRAN incoming feasibility note for a new
submission. The second NOTE is local only:

`checking for future file timestamps ... NOTE`

`unable to verify current time`

No future file timestamps were reported.

## Downstream and additional checks

- Test suite run with `testthat::test_dir()`
- Source archive built with `R CMD build` from a clean copy of the package
- `R CMD check --as-cran` run on the source archive
- URL check via `urlchecker::url_check()`
- Experimental leave-one-out `Di` diagnostics are excluded from this CRAN
  release branch and remain only in the GitHub development branch until the
  methodology is ready for release.

## Resubmission

Initial submission.
