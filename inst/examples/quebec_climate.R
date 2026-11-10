# Real Quebec climate fields. Run after installing GLBFP from develop.
# Rscript quebec_climate.R [output-directory], or source() for a temp directory.
library(GLBFP)
library(ggplot2)
library(plotly)
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1] else file.path(tempdir(), "GLBFP-quebec-climate")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)
data_path <- system.file("extdata", "cangrd_quebec_1991_2020.csv", package = "GLBFP")
if (!nzchar(data_path)) stop("Install the develop version with the bundled CANGRD data.")
climate <- read.csv(data_path)
stopifnot(nrow(climate) == 2760L, length(unique(climate$year)) == 30L,
          nrow(unique(climate[c("i", "j")])) == 92L,
          !anyDuplicated(climate[c("i", "j", "year")]))
values <- c("tmin_anomaly_c", "tmax_anomaly_c")
make_field <- function(table, with_time) {
  glbfp_field(as.matrix(table[values]),
    indices = as.matrix(table[if (with_time) c("i", "j", "year") else c("i", "j")]),
    units = c("deg C", "deg C"), dependence = "dependent")
}
field2 <- make_field(climate[climate$year == 2020L, ], FALSE)
field3 <- make_field(climate, TRUE)
stopifnot(field2$N == 2L, field3$N == 3L, field2$d == 2L, field3$d == 2L)

# Only the evaluation grid is fine. Counting cells remain b/m = 0.2 deg C.
axis <- seq(-3, 5, length.out = 481L)
points <- as.matrix(expand.grid(Tmin = axis, Tmax = axis))
fits <- lapply(list(field2, field3), function(field) {
  GLBFP_estimate(field, b = c(.6, .6), m = c(3L, 3L),
    min_vals = c(-8, -8), max_vals = c(8, 8), grid_points = points,
    engine = "cpp_prefix", inference = "none")
})
caption <- "Data source: Environment and Climate Change Canada, CANGRD v2023.\nHomogenized station observations interpolated by ECCC; no additional spatial interpolation."
map_years <- c(1991L, 2005L, 2020L)
map_limits <- range(climate$tmin_anomaly_c[climate$year %in% map_years])
climate$x_km <- (climate$i - min(climate$i)) * 50
climate$y_km <- (climate$j - min(climate$j)) * 50
map_plot <- function(years, N) {
  ggplot(climate[climate$year %in% years, ], aes(x_km, y_km, fill = tmin_anomaly_c)) +
    geom_tile(width = 50, height = 50) + coord_equal() + facet_wrap(~year, nrow = 1) +
    scale_fill_viridis_c(name = "Tmin anomaly (deg C)", limits = map_limits) +
    labs(title = paste0("Quebec climate field: N = ", N, ", d = 2"),
      subtitle = "Native 50-km grid cells; anomalies relative to 1961-1990",
      x = "Along native i-axis (km)", y = "Along native j-axis (km)", caption = caption) +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom", panel.grid = element_blank(),
      plot.caption = element_text(size = 8.5))
}
save_density_plot(map_plot(2020L, 2L), file.path(output_dir, "N2-climate-map.png"), width = 8, height = 6)
save_density_plot(map_plot(map_years, 3L), file.path(output_dir, "N3-climate-maps.png"), width = 12, height = 5.5)
for (j in seq_along(fits)) {
  N <- j + 1L
  fit <- fits[[j]]
  period <- if (N == 2L) "2020: 92 sites" else "1991-2020: 2,760 site-years pooled"
  labels <- c("Annual Tmin anomaly", "Annual Tmax anomaly")
  plot <- plot_density(fit, "heatmap", labels = labels, units = c("deg C", "deg C"),
    title = paste0("GLBFP density: N = ", N, ", d = 2")) +
    scale_fill_viridis_c(name = "Density (1/deg C squared)") +
    coord_fixed(xlim = c(-3, 5), ylim = c(-3, 5), expand = FALSE) +
    labs(subtitle = period, caption = caption) + theme(plot.caption = element_text(size = 8.5))
  save_density_plot(plot, file.path(output_dir, paste0("N", N, "-climate-density.png")), width = 8, height = 6.5)
  contour <- plot_density(fit, "contour", labels = labels, title = paste("GLBFP contours:", period)) +
    coord_equal() + labs(caption = caption) + theme(plot.caption = element_text(size = 8.5))
  save_density_plot(contour, file.path(output_dir, paste0("N", N, "-climate-contours.png")), width = 8, height = 6.5)
  surface <- plot_density(fit, "surface", labels = labels,
    title = paste0("GLBFP bivariate density<br>", period))
  save_density_plot(surface, file.path(output_dir, paste0("N", N, "-climate-surface.html")))
  write.csv(data.frame(fit$grid, density = fit$densities),
    file.path(output_dir, paste0("N", N, "-climate-density.csv")), row.names = FALSE)
}
write.csv(data.frame(N = c(2L, 3L), d = 2L, observations = c(92L, 2760L),
  b = .6, m = 3L, counting_step = .2, display_step = 8/480,
  anchor = -8, inference = "none"), file.path(output_dir, "example-settings.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Figures and numerical grids written to ", output_dir)
