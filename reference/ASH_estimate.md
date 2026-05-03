# ASH density estimation on a grid

Computes ASH density estimates on a regular or user-supplied grid.

## Usage

``` r
ASH_estimate(
  data,
  b = compute_bi_optim(data, m = rep(1, ncol(data))),
  m = rep(1, ncol(data)),
  grid_size = 20,
  grid_points = NULL,
  min_vals = apply(data, 2, min),
  max_vals = apply(data, 2, max)
)

# S3 method for class 'ASH_estimate'
print(x, ...)

# S3 method for class 'ASH_estimate'
plot(x, contour = FALSE, ...)
```

## Arguments

- data:

  Numeric matrix or data frame of observations (`n x d`).

- b:

  Positive numeric vector of bandwidths (length `d`).

- m:

  Positive integer vector of shifts (length `d`).

- grid_size:

  Integer number of grid points per dimension when `grid_points = NULL`.

- grid_points:

  Optional matrix/data frame of explicit evaluation points.

- min_vals:

  Numeric vector of lower grid bounds (length `d`).

- max_vals:

  Numeric vector of upper grid bounds (length `d`).

- x:

  Object from `ASH_estimate()` to print.

- ...:

  Additional arguments (unused).

- contour:

  If `TRUE`, draw a contour-like 2D representation for 2D data.

## Value

A list with class `c("glbfp_grid", "ASH_estimate")` containing grid
coordinates, densities, and grid metadata.

## Details

When `grid_points` is `NULL`, a regular grid is constructed from
`min_vals` to `max_vals`. Custom grids may be irregular; in that case
plotting uses point or scatter representations instead of a surface.

## Methods (by generic)

- `print(ASH_estimate)`: Print method for object of class
  `"ASH_estimate"`.

- `plot(ASH_estimate)`: Plot method for object of class
  `"ASH_estimate"`.

## See also

[`ASH()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/ASH.md),
[`LBFP_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/LBFP_estimate.md),
[`GLBFP_estimate()`](https://aureliennicosiaulaval.github.io/GLBFP/reference/GLBFP_estimate.md)

## Examples

``` r
b <- c(0.5, 0.5)
# Use a small, representative subset so examples remain fast in checks.
sample_data <- as.matrix(ashua[seq_len(120), -3])
out <- ASH_estimate(sample_data, b = b, m = c(1, 1), grid_size = 10)
out
#> ASH Density Estimation on Grid:
#> Grid points: 100 
#> Dimensions: 2 
#> Bandwidths (b): 0.5, 0.5 
#> Shifts (m): 1, 1 
```
