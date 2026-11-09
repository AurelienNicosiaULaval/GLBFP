#' Visualize estimated densities in one to three observation dimensions
#'
#' A surface for `d=2` uses height for density. For `d=3`, all three spatial
#' plot axes represent observation components; density is shown by slices or
#' isosurfaces. Neither plot dimension refers to the index dimension `N`.
#' @param x A grid fit from [GLBFP_estimate()] or its ASH/LBFP counterparts.
#' @param type `"line"`, `"heatmap"`, `"contour"`, `"surface"`, `"slices"`,
#'   or `"isosurface"`. Default depends on the observation dimension.
#' @param labels Optional vector of observation axis labels (length `d`).
#' @param units Optional vector of units (length `d`); defaults to field units.
#'   Density units are the reciprocal of their product, not the original units.
#' @param slice_dim Observation dimension held fixed for slices (default 3).
#' @param slices Values on that grid axis. Default uses up to three evenly
#'   spaced interior grid planes. A slice is a joint-density section, not a
#'   marginal or conditional density and is not renormalized.
#' @param isomin,isomax Optional lower and upper isosurface density levels.
#' @param surface_count Number of equally spaced isosurfaces (default 3).
#' @param title Optional title.
#' @return A `ggplot` for lines, heatmaps, contours and slices, or a `plotly`
#'   widget for surfaces/isosurfaces. Surfaces and volume plots require distinct
#'   complete Cartesian grids and at most 250,000 points. Unstructured 2D
#'   grids use points for legacy surface/contour requests; heatmaps reject them.
#' @examples
#' set.seed(20)
#' x <- cbind(a = rnorm(100), b = rnorm(100))
#' fit <- GLBFP_estimate(x, b = c(1, 1), grid_size = 12)
#' plot_density(fit, "heatmap", units = c("cm", "s"))
#' plot_density(fit, "contour")
#' @export
plot_density <- function(x, type = NULL, labels = NULL, units = x$units,
                          slice_dim = 3L, slices = NULL, isomin = NULL, isomax = NULL,
                          surface_count = 3L, title = NULL) {
  if (!inherits(x, "glbfp_grid")) stop("`x` must be a grid density fit.", call. = FALSE)
  d <- ncol(x$grid)
  if (d > 3L) stop("Plotting requires d <= 3; fit explicit marginal distributions for larger d.", call. = FALSE)
  if (is.null(type)) type <- c("line", "heatmap", "slices")[d]
  type <- match.arg(type, c("line", "heatmap", "contour", "surface", "slices", "isosurface"))
  valid <- list("line", c("heatmap", "contour", "surface"), c("slices", "isosurface"))
  if (!type %in% valid[[d]]) stop("Plot type does not match the observation dimension d.", call. = FALSE)
  if (is.null(labels)) labels <- x$col_names
  if (is.null(labels)) labels <- paste0("X", seq_len(d))
  if (!is.character(labels) || length(labels) != d || anyNA(labels)) stop("`labels` must have length d.", call. = FALSE)
  density_label <- "Density"
  if (!is.null(units)) {
    if (!is.character(units) || length(units) != d || anyNA(units)) stop("`units` must have length d.", call. = FALSE)
    labels <- paste0(labels, " [", units, "]")
    unit_terms <- vapply(unique(units), function(unit) {
      k <- sum(units == unit)
      if (k == 1L) unit else paste0(unit, "^", k)
    }, character(1))
    density_label <- paste0("Density [1/(", paste(unit_terms, collapse = " "), ")]")
  }
  if (is.null(title)) title <- paste(x$method, "density")
  if (nrow(x$grid) > 250000L) stop("Plot exceeds 250,000 points; evaluate a smaller visualization grid.", call. = FALSE)
  if (d == 1L) {
    df <- data.frame(x = x$grid[, 1], density = x$densities)
    df <- df[order(df$x), ]
    return(ggplot2::ggplot(df, ggplot2::aes(x, density)) + ggplot2::geom_line(linewidth = .7, color = "#176B87") +
      ggplot2::labs(x = labels[1], y = density_label, title = title) + ggplot2::theme_minimal(base_size = 12))
  }
  regular <- !anyDuplicated(as.data.frame(x$grid)) &&
    prod(vapply(seq_len(d), function(j) length(unique(x$grid[, j])), integer(1))) == nrow(x$grid)
  if (!regular) {
    if (d == 3L || type == "heatmap") stop("This plot requires a complete Cartesian grid without duplicate points.", call. = FALSE)
    df <- data.frame(x = x$grid[, 1], y = x$grid[, 2], density = x$densities)
    if (type == "contour") return(ggplot2::ggplot(df, ggplot2::aes(x, y, color = density)) +
      ggplot2::geom_point() + ggplot2::scale_color_viridis_c(name = density_label) +
      ggplot2::labs(x = labels[1], y = labels[2], title = title) + ggplot2::theme_minimal(base_size = 12))
    return(plotly::plot_ly(df, x = ~x, y = ~y, z = ~density, color = ~density,
      type = "scatter3d", mode = "markers") |>
      plotly::layout(title = list(text = title, x = .04, y = .98),
        margin = list(t = 65, b = 45, l = 45, r = 45), scene = list(xaxis = list(title = labels[1]),
        yaxis = list(title = labels[2]), zaxis = list(title = density_label))))
  }
  if (d == 2L) {
    df <- data.frame(x = x$grid[, 1], y = x$grid[, 2], density = x$densities)
    if (type == "surface") {
      xs <- sort(unique(df$x)); ys <- sort(unique(df$y))
      # Plotly expects rows in y, columns in x. Sorting makes row permutations
      # harmless and works on rectangular grids with unequal axis lengths.
      z <- matrix(df$density[order(df$x, df$y)], nrow = length(ys), ncol = length(xs))
      return(plotly::plot_ly(x = xs, y = ys, z = z, type = "surface", colorscale = "Viridis",
         colorbar = list(title = density_label)) |>
        plotly::layout(title = list(text = title, x = .04, y = .98),
        margin = list(t = 65, b = 45, l = 45, r = 45), scene = list(xaxis = list(title = labels[1]),
           yaxis = list(title = labels[2]), zaxis = list(title = density_label))))
    }
    if (type == "contour") return(ggplot2::ggplot(df, ggplot2::aes(x, y, z = density)) +
      ggplot2::geom_contour_filled(bins = 9) + ggplot2::scale_fill_viridis_d(name = density_label) +
      ggplot2::labs(x = labels[1], y = labels[2], title = title) + ggplot2::theme_minimal(base_size = 12))
    return(glbfp_heatmap(df, labels[1:2], density_label, title))
  }
  if (type == "isosurface") {
    if (is.null(isomin)) isomin <- max(x$densities) * .15
    if (is.null(isomax)) isomax <- max(x$densities) * .8
    if (length(isomin) != 1L || length(isomax) != 1L || !is.numeric(isomin) ||
        !is.numeric(isomax) || !is.finite(isomin) || !is.finite(isomax) ||
        isomin < 0 || isomin >= isomax || isomax > max(x$densities)) {
      stop("Isosurface levels must satisfy 0 <= isomin < isomax <= maximum density.", call. = FALSE)
    }
    surface_count <- glbfp_validate_vector(surface_count, 1, "surface_count", positive = TRUE, integerish = TRUE)
    return(plotly::plot_ly(x = x$grid[, 1], y = x$grid[, 2], z = x$grid[, 3], value = x$densities,
      type = "isosurface", isomin = isomin, isomax = isomax,
      surface = list(count = surface_count), opacity = .45,
      caps = list(x = list(show = FALSE), y = list(show = FALSE), z = list(show = FALSE)),
      colorscale = "Viridis", colorbar = list(title = density_label)) |>
      plotly::layout(title = list(text = title, x = .04, y = .98),
        margin = list(t = 65, b = 45, l = 45, r = 45), scene = list(xaxis = list(title = labels[1]),
        yaxis = list(title = labels[2]), zaxis = list(title = labels[3]))))
  }
  slice_dim <- glbfp_validate_vector(slice_dim, 1, "slice_dim", positive = TRUE, integerish = TRUE)
  if (slice_dim > 3L) stop("`slice_dim` must be 1, 2 or 3.", call. = FALSE)
  axis <- sort(unique(x$grid[, slice_dim]))
  if (is.null(slices)) slices <- unique(axis[pmax(1L, round(seq(1, length(axis), length.out = 5)[2:4]))])
  if (!is.numeric(slices) || !length(slices) || any(!is.finite(slices)) || !all(slices %in% axis)) {
    stop("`slices` must be values on the selected grid axis.", call. = FALSE)
  }
  remaining <- setdiff(1:3, slice_dim)
  keep <- x$grid[, slice_dim] %in% slices
  df <- data.frame(x = x$grid[keep, remaining[1]], y = x$grid[keep, remaining[2]],
    density = x$densities[keep], plane = paste0(labels[slice_dim], " = ",
       format(signif(x$grid[keep, slice_dim], 4), trim = TRUE)))
  glbfp_heatmap(df, labels[remaining], density_label, title) + ggplot2::facet_wrap(~plane)
}

glbfp_heatmap <- function(df, labels, density_label, title) {
  # Midpoint edges also handle nonuniform Cartesian evaluation grids.
  for (variable in c("x", "y")) {
    axis <- sort(unique(df[[variable]]))
    if (length(axis) < 2L) stop("Heatmap axes require at least two grid coordinates.", call. = FALSE)
    edges <- c(axis[1] - diff(axis)[1] / 2, (utils::head(axis, -1) + utils::tail(axis, -1)) / 2,
               utils::tail(axis, 1) + utils::tail(diff(axis), 1) / 2)
    at <- match(df[[variable]], axis)
    df[[paste0(variable, "min")]] <- edges[at]
    df[[paste0(variable, "max")]] <- edges[at + 1L]
  }
  ggplot2::ggplot(df, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = density)) +
    ggplot2::geom_rect() + ggplot2::scale_fill_viridis_c(name = density_label) +
    ggplot2::labs(x = labels[1], y = labels[2], title = title) + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
      barwidth = grid::unit(7, "cm"), barheight = grid::unit(.35, "cm"))) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")
}

#' Export a density visualization
#' @param plot A ggplot or plotly object from [plot_density()].
#' @param file Output path: PNG, PDF or SVG for ggplot, HTML for plotly.
#' @param width,height Figure dimensions in inches for static output.
#' @param dpi Raster resolution.
#' @return The output path invisibly. HTML is self-contained.
#' @export
save_density_plot <- function(plot, file, width = 8, height = 6, dpi = 150) {
  extension <- tolower(tools::file_ext(file))
  if (inherits(plot, "ggplot") && extension %in% c("png", "pdf", "svg")) {
    ggplot2::ggsave(file, plot = plot, width = width, height = height, dpi = dpi)
  } else if (inherits(plot, "plotly") && extension == "html") {
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("HTML export requires htmlwidgets.", call. = FALSE)
    htmlwidgets::saveWidget(plot, file, selfcontained = TRUE)
  } else stop("Use PNG/PDF/SVG for static plots or HTML for interactive plots.", call. = FALSE)
  invisible(file)
}
