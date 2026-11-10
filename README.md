<img src="man/figures/logo.svg" align="right" width="160" alt="GLBFP hex logo" />

# GLBFP

[![R-CMD-check](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml/badge.svg?branch=develop)](https://github.com/AurelienNicosiaULaval/GLBFP/actions/workflows/R-CMD-check.yaml?query=branch%3Adevelop)
[![CRAN status](https://www.r-pkg.org/badges/version/GLBFP)](https://CRAN.R-project.org/package=GLBFP)
[![License](https://img.shields.io/cran/l/GLBFP)](https://CRAN.R-project.org/package=GLBFP)

`GLBFP` implements Averaged Shifted Histograms (ASH), Linear Blend Frequency
Polygons (LBFP), and General Linear Blend Frequency Polygons (GLBFP) for
nonparametric density estimation. It provides pointwise and grid evaluation,
reusable native sparse counts, indexed fields, fixed-grid leave-one-out
self-support scores, and density plots for one to three observation variables.

## Version and installation

| Source | Version | Status |
| --- | --- | --- |
| CRAN | 0.5.2 | Stable release; no `D_i` diagnostics |
| GitHub `main` | 0.5.2.9000 | Earlier development version |
| GitHub `develop` | 0.5.2.9001 | Development version for scientific review |

To examine the corrections, indexed fields and real Quebec examples described
here, install `develop` explicitly. Source installation requires a C++ compiler.

```r
install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/GLBFP", ref = "develop")
library(GLBFP)
packageVersion("GLBFP")
```

The stable CRAN release remains available with `install.packages("GLBFP")`.
The [scientific review guide](review/README.md) records the scope, examples,
validation evidence and limitations of 0.5.2.9001. This version is not a CRAN
release and does not replace the submitted article's replication software.

## Real Quebec data: N=2 and N=3, with d=2

The bundled CANGRD subset contains 92 Quebec grid sites for each year from
1991 to 2020. ECCC homogenized station observations and interpolated them onto
its native 50-km lattice. No further spatial interpolation is applied here.
The two measured components are anomalies of annual mean daily minimum and
maximum temperatures relative to 1961-1990, in degrees Celsius.

- `N=2`: the spatial indices `(i,j)` for 2020, with 92 observations.
- `N=3`: the spatial and year indices `(i,j,year)`, with 2,760 observations.
- `d=2` in both cases: the two temperature-anomaly variables.

```r
library(GLBFP)
climate <- read.csv(system.file(
  "extdata", "cangrd_quebec_1991_2020.csv", package = "GLBFP"
))
field <- glbfp_field(
  values = as.matrix(climate[c("tmin_anomaly_c", "tmax_anomaly_c")]),
  indices = as.matrix(climate[c("i", "j", "year")]),
  units = c("deg C", "deg C"), dependence = "dependent"
)
field
points <- as.matrix(expand.grid(
  tmin = seq(-3, 5, length.out = 201),
  tmax = seq(-3, 5, length.out = 201)
))
fit <- GLBFP_estimate(
  field, b = c(.6, .6), m = c(3, 3),
  min_vals = c(-8, -8), max_vals = c(8, 8), grid_points = points
)
plot_density(fit, "heatmap", labels = c("Tmin anomaly", "Tmax anomaly"))
plot_density(fit, "contour", labels = c("Tmin anomaly", "Tmax anomaly"))
plot_density(fit, "surface", labels = c("Tmin anomaly", "Tmax anomaly"))
```

The density describes the observed anomaly pairs, with equal weights per
selected site-year. It does not estimate a density over spatial indices.
A 3D surface uses density as height; the observation dimension remains `d=2`.
Native geographic grid cells remain visible in the maps. The smoothness of a
density display is controlled separately by its evaluation grid.

Run the installed example to export maps, heatmaps, contours and interactive
surfaces without downloading any additional data:

```r
source(system.file("examples", "quebec_climate.R", package = "GLBFP"))
```

The [climate vignette](vignettes/quebec-climate-fields.Rmd) explains the data
structure, grid settings, unequal roles of `N` and `d`, and interpretation.
Source hashes, coordinates, extraction details and data-licence attribution
are included in `inst/extdata/cangrd_quebec_PROVENANCE.txt` and
`inst/extdata/cangrd_quebec_SOURCE.json`. A separate CDEM elevation example
near Quebec City is also included. These applications are descriptive.

## Grid conventions and corrected behaviour

For bandwidth `b`, integer `m`, and anchor `a=min_vals`, fine counting cells
are half-open intervals `[a+k*b/m, a+(k+1)*b/m)`. Interpolation uses cell
centres. The lattice is unbounded and unoccupied cells have zero counts.
`max_vals` limits a plotting grid; it does not clamp or renormalize density
support. The previous tight-bound extrapolation is corrected.

`compute_bi_optim()` retains its historical name and implements a
Gaussian-reference L1 bandwidth heuristic from Carbon's working note.
Two transcription errors are corrected and the requested `m` is used.
Automatically selected widths and affected boundary values intentionally
change. Historical calls, uppercase functions, lowercase aliases and legacy
classes remain available.

`glbfp_field()` accepts explicit integer sites or scalar/multicomponent
arrays, stores `N` separately from `d`, rejects duplicate sites, records omitted
sites and accepts incomplete grids without imputation. A field object does
not establish stationarity or mixing assumptions.

## Computation and fixed-grid diagnostics

```r
model <- prepare_glbfp(field, b = c(.6, .6), m = c(3, 3),
                       min_vals = c(-8, -8), max_vals = c(8, 8))
predict(model, matrix(c(0, 0), nrow = 1))
scores <- compute_Di(model)
summary(scores)
```

`prepare_glbfp()` reuses the same counts for predictions and all `D_i` scores.
`compute_Di()` is the only public name for this score. The grid and bandwidth
remain fixed when each observation is deleted. The score measures local
self-support; it is not an outlier probability or a calibrated test.

`engine="auto"` selects native direct enumeration for stencils of at most
64 cells and native prefix traversal otherwise. `R`, `cpp_direct` and
`cpp_prefix` can be selected explicitly. Gains depend on sample size, dimension,
occupancy, `m`, query count and reuse of preparation. Prefix traversal retains
a combinatorial worst case. The [review guide](review/README.md) reports
favourable and unfavourable measured regimes. MIT-licensed native adaptations
and their notices are retained in `inst/COPYRIGHTS`; the package is GPL (>= 3).

## Inference and scope

Uncertainty defaults to `inference="none"`, with `sd` and `IC` set to `NA`.
Opt-in `inference="iid"` requires independent observations, explicit fixed
bandwidths and anchor, and evaluation points fixed independently of the data.
Its variance is the sample variance of observation weights divided by `n`.
Normal intervals target the expected smoothed estimator, with smoothing bias
remaining. There is no implemented dependent-field confidence interval.

The single-grid LBFP diagnostic studied by Nicosia, Duchesne and Carbon
(2026, [arXiv:2607.14335](https://arxiv.org/abs/2607.14335)) has a separate,
frozen `LOOLBFP` replication package. Use that environment to reproduce the
submitted JCGS results. Equivalent fixed grids give agreement to floating-point
precision; default anchors differ. The package's ASH/GLBFP diagnostic extensions
and arbitrary index representation do not inherit all results of that article.

## Visualization and documentation

`plot_density()` supports lines for `d=1`, heatmaps/contours/surfaces for `d=2`,
and joint-density slices/isosurfaces for `d=3`. `save_density_plot()` exports
PNG, PDF, SVG, or self-contained interactive HTML. Simulated examples with
known marginal laws are retained for numerical validation in
`inst/examples/indexed_fields.R`.

The package includes function manuals and executable vignettes on estimator
choice, indexed fields, computation, diagnostics and validation. Built vignettes
are included in source archives attached to the review release.

## References

- Carbon, M. and Duchesne, T. (2024). Multivariate frequency polygon for
  stationary random fields. Annals of the Institute of Statistical Mathematics,
  76, 263-287. <https://doi.org/10.1007/s10463-023-00883-5>.
- Carbon, M. and Duchesne, T. (2026; online 2025). Asymptotic normality of
  multivariate frequency polygons for stationary random fields. Annals of the
  Institute of Statistical Mathematics, 78, 297-326.
  <https://doi.org/10.1007/s10463-025-00952-x>.
- Carbon, M. (2024). Le fameux GLBFP. Internal working note, Universite Laval,
  equations (2)-(7) and (22)-(36). This note specifies the package construction
  and Gaussian-reference L1 heuristic; general GLBFP normality is left unproved.
- Scott, D. W. (1992). Multivariate Density Estimation: Theory, Practice, and
  Visualization. Wiley. <https://doi.org/10.1002/9780470316849>.
- Terrell, G. R. and Scott, D. W. (1985). Oversmoothed Nonparametric Density
  Estimates. Journal of the American Statistical Association, 80, 209-214.
  <https://doi.org/10.1080/01621459.1985.10477163>.

Use `citation("GLBFP")` for software citation metadata. Report the version and
commit used. The software concept DOI is <https://doi.org/10.5281/zenodo.17945962>;
it does not establish that every development commit has been archived.

## Development checks

```r
devtools::test()
devtools::check()
```

Benchmarks are stored in `benchmarks/` and are not run during `R CMD check`.
Please report issues at <https://github.com/AurelienNicosiaULaval/GLBFP/issues>.
