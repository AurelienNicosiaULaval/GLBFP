# Optional data acquisition, run from the GLBFP source-package root. Downloads a small
# window via GDAL HTTP range reads. Never crops a projected raster in degrees.
library(terra)
stopifnot(file.exists("DESCRIPTION"), read.dcf("DESCRIPTION")[1, "Package"] == "GLBFP")
url <- "https://datacube-prod-data-public.s3.ca-central-1.amazonaws.com/store/elevation/cdem-cdsm/cdem/cdem-canada-dem.tif"
dem <- rast(paste0("/vsicurl/", url))
stopifnot(grepl("3979", crs(dem)), all(res(dem) == 16))
center <- project(vect(data.frame(lon = -71.305, lat = 46.805),
                       geom = c("lon", "lat"), crs = "EPSG:4326"), crs(dem))
xy <- crds(center)[1, ]
window <- crop(dem, ext(xy[1] - 256, xy[1] + 256, xy[2] - 256, xy[2] + 256), snap = "near")
table <- as.data.frame(window, xy = TRUE, cells = TRUE, na.rm = FALSE)
names(table)[4] <- "elevation_m"
rc <- rowColFromCell(window, table$cell)
table$i <- rc[, 2]
table$j <- nrow(window) + 1L - rc[, 1]
stopifnot(nrow(table) >= 500, all(is.finite(table$elevation_m)),
          !anyDuplicated(table[c("i", "j")]), terra::is.lonlat(window) == FALSE)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
write.csv(table[c("i", "j", "x", "y", "elevation_m")],
          "inst/extdata/cdem_quebec.csv", row.names = FALSE)
writeRaster(window, "data-raw/cdem-quebec-3979.tif", overwrite = TRUE)
writeLines(c(paste("Source:", url), paste("Retrieved:", Sys.Date()),
  "Product: CDEM; publisher: Natural Resources Canada",
  "Catalogue: https://open.canada.ca/data/en/dataset/7f245e4d-76c2-4caa-951a-45d1d2051333",
  "Licence: Open Government Licence - Canada",
  "https://open.canada.ca/en/open-government-licence-canada",
  "Contains information licensed under the Open Government Licence - Canada.",
  "Small window near Quebec City; center requested in EPSG:4326 (-71.305, 46.805).",
  "Coordinates transformed to the source CRS before cropping. No resampling.",
  "CRS: EPSG:3979; source COG grid resolution: 16 m; elevation units: m.",
  paste("Rows:", nrow(window), "Columns:", ncol(window)),
  paste("Source extent:", paste(as.vector(ext(window)), collapse = ", ")),
  "i increases eastward, j northward; x/y are projected cell-center coordinates.",
  "This is a descriptive topographic distribution. Stationarity/mixing are not asserted.",
  "CDEM source acquisition dates vary by location; extraction date is not a measurement date.",
  "No endorsement by Natural Resources Canada is implied."),
  "inst/extdata/cdem_quebec_PROVENANCE.txt")
writeLines(capture.output(list(raster = window, gdal = terra::gdal(), session = sessionInfo())),
           "data-raw/cdem-extraction-environment.txt")
