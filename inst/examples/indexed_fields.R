# Execute after installing GLBFP. Optional first argument: output directory.
library(GLBFP)
library(ggplot2)
library(plotly)
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1] else file.path(tempdir(), "GLBFP-indexed-fields")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)
set.seed(20261109)

# Stationary finite-range Gaussian moving average on Z^N. Innovations outside
# the observation box are simulated as well, so boundary sites have the same
# marginal distribution. Different observation components use independent
# innovations before applying the specified covariance matrix.
simulate_ma_field <- function(shape, coefficients, Sigma) {
  N <- length(shape); d <- ncol(Sigma)
  stopifnot(length(coefficients) == N, nrow(Sigma) == d)
  sites <- as.matrix(expand.grid(lapply(shape, seq_len)))
  values <- matrix(0, nrow(sites), d)
  for (component in seq_len(d)) {
    innovations <- array(rnorm(prod(shape + 1L)), shape + 1L)
    value <- innovations[sites + 1L]
    for (j in seq_len(N)) {
      neighbor <- sites + 1L
      neighbor[, j] <- neighbor[, j] - 1L
      value <- value + coefficients[j] * innovations[neighbor]
    }
    values[, component] <- value / sqrt(1 + sum(coefficients^2))
  }
  values <- values %*% chol(Sigma)
  colnames(values) <- paste0("X", seq_len(d))
  glbfp_field(values, sites, units = rep("standard unit", d), dependence = "dependent")
}

normal_density <- function(points, Sigma) {
  exp(-rowSums((points %*% solve(Sigma)) * points) / 2) /
    sqrt((2*pi)^ncol(points) * det(Sigma))
}

# N=1 remains the existing observations-in-rows workflow.
field1 <- simulate_ma_field(200L, .5, matrix(1))
grid1 <- GLBFP_estimate(field1, b = .65, m = 3, min_vals = -5, max_vals = 5, grid_size = 150)
p1 <- plot_density(grid1) + geom_line(data = data.frame(x = grid1$grid[, 1], density = dnorm(grid1$grid[, 1])),
                                      aes(x, density), color = "#BA4A00", linetype = 2) +
  labs(subtitle = "Solid: GLBFP estimate; dashed: known N(0,1) marginal")
save_density_plot(p1, file.path(output_dir, "N1-density.png"))

# N=2, d=2. The joint marginal law is exactly Gaussian with this covariance.
Sigma2 <- matrix(c(1, .45, .45, 1), 2, 2)
field2 <- simulate_ma_field(c(35, 30), c(.6, .4), Sigma2)
map2 <- data.frame(field2$indices, X1 = field2$values[, 1])
names(map2)[1:2] <- c("i", "j")
pmap <- ggplot(map2, aes(i, j, fill = X1)) + geom_raster() + coord_equal() +
  scale_fill_viridis_c(name = "X1 [standard unit]") + theme_minimal(base_size = 12) +
  labs(title = "Simulated field: N = 2, d = 2", x = "Index i", y = "Index j")
save_density_plot(pmap, file.path(output_dir, "N2-field.png"))
fit2 <- GLBFP_estimate(field2, b = c(.65, .65), m = c(3, 3), grid_size = 55,
                       min_vals = c(-5, -5), max_vals = c(5, 5))
truth2 <- normal_density(fit2$grid, Sigma2)
comparison <- rbind(data.frame(fit2$grid, density = fit2$densities, source = "GLBFP estimate"),
                    data.frame(fit2$grid, density = truth2, source = "Known Gaussian marginal"))
pcompare <- ggplot(comparison, aes(X1, X2, fill = density)) + geom_raster() + coord_equal() +
  facet_wrap(~source) + scale_fill_viridis_c(name = "Joint density [1/unit²]") + theme_minimal(base_size = 12) +
  labs(title = "Estimate and known marginal, N = 2, d = 2", x = "X1 [standard unit]", y = "X2 [standard unit]") +
  theme(legend.position = "bottom")
save_density_plot(pcompare, file.path(output_dir, "N2-truth-comparison.png"), width = 10, height = 5)
save_density_plot(plot_density(fit2, "contour"), file.path(output_dir, "N2-contours.png"))
save_density_plot(plot_density(fit2, "heatmap"), file.path(output_dir, "N2-heatmap.png"))
save_density_plot(plot_density(fit2, "surface"), file.path(output_dir, "N2-surface.html"))
write.csv(data.frame(fit2$grid, estimate = fit2$densities, truth = truth2),
          file.path(output_dir, "N2-density-grid.csv"), row.names = FALSE)
write.csv(data.frame(field2$indices, field2$values), file.path(output_dir, "N2-observations.csv"), row.names = FALSE)

# N=3, d=3. The isosurface axes are X1, X2, X3, not the three site indices.
Sigma3 <- diag(c(1, 1.3, .7))
field3 <- simulate_ma_field(c(14, 12, 10), c(.5, .3, .2), Sigma3)
fit3 <- GLBFP_estimate(field3, b = c(.8, .9, .7), m = c(2, 2, 2), grid_size = 19,
                       min_vals = c(-4, -4, -4), max_vals = c(4, 4, 4))
truth3 <- normal_density(fit3$grid, Sigma3)
save_density_plot(plot_density(fit3, "slices"), file.path(output_dir, "N3-density-slices.png"), width = 11, height = 4.5)
save_density_plot(plot_density(fit3, "isosurface"), file.path(output_dir, "N3-isosurfaces.html"))
plane <- field3$indices[, 3] == 5L
df <- data.frame(i = field3$indices[plane, 1], j = field3$indices[plane, 2], X1 = field3$values[plane, 1])
spatial3 <- ggplot(df, aes(i, j, fill = X1)) + geom_raster() + coord_equal() +
  scale_fill_viridis_c(name = "X1 [standard unit]") + theme_minimal(base_size = 12) +
  labs(title = "N = 3 field: index plane k = 5", x = "Index i", y = "Index j")
save_density_plot(spatial3, file.path(output_dir, "N3-field-plane.png"))
write.csv(data.frame(fit3$grid, estimate = fit3$densities, truth = truth3),
          file.path(output_dir, "N3-density-grid.csv"), row.names = FALSE)
write.csv(data.frame(field3$indices, field3$values), file.path(output_dir, "N3-observations.csv"), row.names = FALSE)

# A real N=2 lattice: a 16 m CDEM window near Quebec City (EPSG:3979).
# Attribution, source URL, extraction and limits are shipped beside the CSV.
file <- system.file("extdata", "cdem_quebec.csv", package = "GLBFP")
dem <- read.csv(file)
real_field <- glbfp_field(matrix(dem$elevation_m, ncol = 1, dimnames = list(NULL, "Elevation")),
                          as.matrix(dem[c("i", "j")]), units = "m", index_units = c("16 m", "16 m"),
                          dependence = "unspecified")
real_fit <- GLBFP_estimate(real_field, b = 3, m = 3, grid_size = 200,
  min_vals = min(dem$elevation_m) - 5, max_vals = max(dem$elevation_m) + 5)
save_density_plot(plot_density(real_fit, title = "CDEM window: descriptive elevation density"),
                   file.path(output_dir, "quebec-elevation-density.png"))
real_map <- ggplot(dem, aes(x, y, fill = elevation_m)) + geom_raster() + coord_equal() +
  scale_fill_viridis_c(name = "Elevation [m]") + theme_minimal(base_size = 12) +
  labs(title = "CDEM near Quebec City", x = "Easting [m, EPSG:3979]", y = "Northing [m, EPSG:3979]") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_density_plot(real_map, file.path(output_dir, "quebec-elevation-map.png"))

diagnostics <- data.frame(example = c("N1_d1", "N2_d2", "N3_d3", "CDEM_N2_d1"),
  n = c(nrow(field1$values), nrow(field2$values), nrow(field3$values), nrow(real_field$values)),
  N = c(1, 2, 3, 2), d = c(1, 2, 3, 1),
  mean_squared_grid_error = c(mean((grid1$densities-dnorm(grid1$grid[, 1]))^2),
    mean((fit2$densities-truth2)^2), mean((fit3$densities-truth3)^2), NA_real_))
write.csv(diagnostics, file.path(output_dir, "examples-summary.csv"), row.names = FALSE)
saveRDS(list(N1 = field1, N2 = field2, N3 = field3, real = real_field,
             fit2 = fit2, fit3 = fit3), file.path(output_dir, "examples.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
