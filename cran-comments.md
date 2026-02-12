## Test environments

- local macOS (aarch64, R 4.5.0)
- GitHub Actions: macOS-latest (release), windows-latest (release), ubuntu-latest (devel/release/oldrel-1)

## R CMD check results

`R CMD check --as-cran` was run on local source tarball.

- 0 ERROR
- 0 WARNING
- 1 NOTE

NOTE is the expected CRAN incoming check note for a new submission.

## Downstream and additional checks

- `devtools::test()`
- `devtools::check(manual = FALSE)`
- URL check via `urlchecker::url_check()`

## Resubmission

Initial submission.
