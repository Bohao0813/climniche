.metric_label <- function(metric) {
  metric <- .metric_key(metric)
  switch(metric,
         climate_change_amount = "Climatic Displacement",
         niche_distance_change = "Niche Distance Shift",
         climate_reconfiguration = "Climatic Reconfiguration",
         change_alignment = "Change alignment",
         niche_boundary_exceedance = "Niche Boundary Exceedance",
         metric)
}

.need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for this plot.", call. = FALSE)
  }
}

.need_raster <- function() {
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("raster is required for raster plotting.", call. = FALSE)
  }
}

.need_terra <- function() {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("terra is required for terra raster plotting.", call. = FALSE)
  }
}

.is_raster <- function(x) {
  methods::is(x, "Raster")
}

.is_spatraster <- function(x) {
  inherits(x, "SpatRaster")
}

.spatial_df <- function(x, value = "value") {
  if (.is_raster(x)) {
    return(.raster_df(x, value = value))
  }
  if (.is_spatraster(x)) {
    return(.terra_df(x, value = value))
  }
  stop("x must be a RasterLayer or terra SpatRaster.", call. = FALSE)
}

.raster_df <- function(x, value = "value") {
  .need_raster()
  if (raster::nlayers(x) != 1L) {
    stop("Only one-layer Raster* objects can be plotted.", call. = FALSE)
  }
  pts <- raster::rasterToPoints(x)
  pts <- as.data.frame(pts)
  names(pts) <- c("x", "y", value)
  pts
}

.terra_df <- function(x, value = "value") {
  .need_terra()
  if (terra::nlyr(x) != 1L) {
    stop("Only one-layer SpatRaster objects can be plotted.", call. = FALSE)
  }
  pts <- terra::as.data.frame(x, xy = TRUE, na.rm = FALSE)
  names(pts)[seq_len(3L)] <- c("x", "y", value)
  pts <- pts[!is.na(pts[[value]]), c("x", "y", value), drop = FALSE]
  pts
}

.spatial_res <- function(x) {
  if (.is_raster(x)) {
    .need_raster()
    return(raster::res(x))
  }
  if (.is_spatraster(x)) {
    .need_terra()
    return(terra::res(x))
  }
  stop("x must be a RasterLayer or terra SpatRaster.", call. = FALSE)
}

.spatial_is_lonlat <- function(x) {
  if (.is_raster(x)) {
    .need_raster()
    return(isTRUE(raster::isLonLat(x)))
  }
  if (.is_spatraster(x)) {
    .need_terra()
    return(isTRUE(terra::is.lonlat(x)))
  }
  FALSE
}

.degree_label <- function(x, positive, negative) {
  value <- abs(x)
  label <- ifelse(
    abs(value - round(value)) < 1e-7,
    as.character(round(value)),
    formatC(value, format = "fg", digits = 4)
  )
  hemisphere <- ifelse(x < 0, negative, ifelse(x > 0, positive, ""))
  paste0(label, "\u00b0", hemisphere)
}

.map_axis_spec <- function(r, degree_labels = c("auto", "none", "hemisphere")) {
  degree_labels <- match.arg(degree_labels)
  use_degrees <- identical(degree_labels, "hemisphere") ||
    (identical(degree_labels, "auto") && .spatial_is_lonlat(r))
  if (!use_degrees) {
    return(list(x = ggplot2::waiver(), y = ggplot2::waiver(),
                xlab = NULL, ylab = NULL))
  }
  list(
    x = function(z) .degree_label(z, "E", "W"),
    y = function(z) .degree_label(z, "N", "S"),
    xlab = "Longitude",
    ylab = "Latitude"
  )
}

.map_extent <- function(dat, extent = NULL) {
  if (is.null(extent)) {
    return(.crop_df(dat))
  }
  extent <- as.numeric(extent)
  if (length(extent) != 4L || any(!is.finite(extent)) ||
      extent[1] >= extent[2] || extent[3] >= extent[4]) {
    stop("extent must be c(xmin, xmax, ymin, ymax).", call. = FALSE)
  }
  list(xlim = extent[1:2], ylim = extent[3:4])
}

.study_region_df <- function(study_region, raster) {
  if (is.null(study_region)) {
    return(NULL)
  }
  if (is.data.frame(study_region) && !inherits(study_region, "sf")) {
    if (!all(c("x", "y") %in% names(study_region))) {
      stop("A study_region data frame must contain x and y columns.",
           call. = FALSE)
    }
    if (!"group" %in% names(study_region)) {
      study_region$group <- 1L
    }
    return(study_region[, c("x", "y", "group"), drop = FALSE])
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("sf is required to draw spatial study_region objects.",
         call. = FALSE)
  }
  region <- try(sf::st_as_sf(study_region), silent = TRUE)
  if (inherits(region, "try-error")) {
    stop("study_region must be an sf, sfc, Spatial, SpatVector, or x-y data frame.",
         call. = FALSE)
  }
  raster_crs <- if (.is_raster(raster)) {
    sf::st_crs(raster::crs(raster))
  } else {
    sf::st_crs(terra::crs(raster, proj = TRUE))
  }
  if (!is.na(sf::st_crs(region)) && !is.na(raster_crs) &&
      sf::st_crs(region) != raster_crs) {
    region <- sf::st_transform(region, raster_crs)
  }
  lines <- suppressWarnings(sf::st_cast(
    sf::st_boundary(sf::st_geometry(region)),
    "LINESTRING"
  ))
  coordinates <- sf::st_coordinates(lines)
  group_columns <- grep("^L[0-9]+$", colnames(coordinates))
  group <- if (length(group_columns)) {
    interaction(as.data.frame(coordinates[, group_columns, drop = FALSE]),
                drop = TRUE)
  } else {
    factor(rep(1L, nrow(coordinates)))
  }
  data.frame(x = coordinates[, "X"], y = coordinates[, "Y"], group = group)
}

.metric_map_colours <- function(metric) {
  switch(
    .metric_key(metric),
    niche_distance_change = c("#386fa4", "#ffffff", "#c7514a"),
    change_alignment = c("#386fa4", "#ffffff", "#c7514a"),
    niche_boundary_exceedance = c("#ffffff", "#f3dfb8", "#d9942f", "#8a3f20"),
    climate_reconfiguration = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45"),
    c("#f7fbff", "#d6e6f2", "#91b9d5", "#27658f")
  )
}

.squish_values <- function(x, range, only.finite = TRUE) {
  if (only.finite) {
    finite <- is.finite(x)
    x[finite] <- pmax(range[1], pmin(range[2], x[finite]))
  } else {
    x <- pmax(range[1], pmin(range[2], x))
  }
  x
}

.mask_to_occupied <- function(x, occupied, occupied_threshold = NULL) {
  if (.is_raster(x)) {
    .need_raster()
    if (!.is_raster(occupied)) {
      stop("occupied must be a RasterLayer when x is a RasterLayer.",
           call. = FALSE)
    }
    if (!raster::compareRaster(x, occupied, stopiffalse = FALSE)) {
      stop("occupied raster must match x geometry.", call. = FALSE)
    }
    mask <- raster::raster(occupied)
    values <- .clean_reference_weights(
      raster::getValues(occupied),
      threshold = occupied_threshold
    )
    raster::values(mask) <- ifelse(values > 0, 1, NA_real_)
    return(raster::mask(x, mask))
  }
  if (.is_spatraster(x)) {
    .need_terra()
    if (!.is_spatraster(occupied)) {
      stop("occupied must be a SpatRaster when x is a SpatRaster.",
           call. = FALSE)
    }
    if (terra::nlyr(occupied) != 1L) {
      stop("occupied must have one layer.", call. = FALSE)
    }
    if (!terra::compareGeom(x, occupied, stopOnError = FALSE)) {
      stop("occupied raster must match x geometry.", call. = FALSE)
    }
    mask <- occupied
    values <- .clean_reference_weights(
      terra::values(mask)[, 1],
      threshold = occupied_threshold
    )
    terra::values(mask) <- ifelse(values > 0, 1, NA_real_)
    return(terra::mask(x, mask))
  }
  stop("x must be a RasterLayer or terra SpatRaster.", call. = FALSE)
}

.occupied_df <- function(occupied, occupied_threshold = NULL) {
  if (is.null(occupied)) {
    return(NULL)
  }
  pts <- .spatial_df(occupied, value = "occupied")
  pts$occupied_weight <- .clean_reference_weights(
    pts$occupied,
    threshold = occupied_threshold
  )
  pts[pts$occupied_weight > 0, , drop = FALSE]
}

.crop_df <- function(dat, pad = 0.02) {
  xr <- range(dat$x, na.rm = TRUE)
  yr <- range(dat$y, na.rm = TRUE)
  dx <- diff(xr) * pad
  dy <- diff(yr) * pad
  list(xlim = c(xr[1] - dx, xr[2] + dx), ylim = c(yr[1] - dy, yr[2] + dy))
}

.plot_title <- function(title, default) {
  if (isFALSE(title)) {
    NULL
  } else if (is.null(title)) {
    default
  } else {
    title
  }
}

.climniche_theme <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0,
                                         colour = "black",
                                         size = base_size + 0.2,
                                         margin = ggplot2::margin(b = 0.8)),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(colour = "black",
                                            size = base_size - 0.8),
      axis.text = ggplot2::element_text(colour = "black",
                                        size = base_size - 1.4),
      axis.title = ggplot2::element_text(colour = "black",
                                         size = base_size - 0.5),
      axis.line = ggplot2::element_line(linewidth = 0.25,
                                        colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25,
                                         colour = "black"),
      legend.position = "right",
      legend.text = ggplot2::element_text(colour = "black",
                                          size = base_size - 1.5),
      legend.title = ggplot2::element_text(colour = "black",
                                           size = base_size - 1.2),
      legend.key.height = grid::unit(4.6, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "black",
                                         size = base_size - 0.4)
    )
}

.map_theme <- function() {
  ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0,
                                         colour = "black", size = 8.7),
      axis.text = ggplot2::element_text(colour = "black", size = 6.4),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.2),
      axis.line = ggplot2::element_line(linewidth = 0.2, colour = "black"),
      legend.position = "right",
      legend.text = ggplot2::element_text(colour = "black", size = 6.8),
      legend.title = ggplot2::element_text(colour = "black", size = 7.2),
      legend.key.width = grid::unit(3.4, "mm"),
      legend.key.height = grid::unit(8, "mm")
    )
}

.climniche_plot_data <- function(x, scope = c("current", "all"),
                                 max_points = 5000, seed = 1) {
  scope <- match.arg(scope)
  max_points <- .check_positive_integer(max_points, "max_points")
  seed <- .check_finite_scalar(seed, "seed")
  if (seed != floor(seed)) {
    stop("seed must be an integer.", call. = FALSE)
  }
  dat <- climniche_table(x, scope = scope)
  if (nrow(dat) > max_points) {
    set.seed(seed)
    dat <- dat[sample(seq_len(nrow(dat)), max_points), , drop = FALSE]
  }
  dat
}

#' Plot a climniche raster metric
#'
#' @param x A `climniche_fit` object with raster outputs, a RasterLayer or a
#'   one-layer terra SpatRaster.
#' @param metric Metric to plot when `x` is a `climniche_fit` object.
#' @param occupied Optional current reference raster to overlay.
#' @param occupied_only If TRUE, mask the plotted raster to current occurrence
#'   cells.
#' @param occupied_threshold Threshold used when `occupied` contains binary or
#'   continuous values. Values above the threshold keep their original value
#'   when used as an overlay or mask.
#' @param title Optional plot title. Use `FALSE` to suppress it.
#' @param midpoint Midpoint for the Niche Distance Shift colour scale.
#' @param limits Optional two-element colour scale limits.
#' @param breaks Optional colour scale breaks.
#' @param colours Optional colour vector replacing the metric palette.
#' @param legend_title Optional colour legend title. Use `FALSE` to suppress
#'   the title.
#' @param legend_position Position passed to the ggplot theme.
#' @param show_legend If FALSE, suppress the colour legend.
#' @param symmetric If TRUE, use limits symmetric around `midpoint`. The
#'   default applies this to Niche Distance Shift and change alignment.
#' @param extent Optional `c(xmin, xmax, ymin, ymax)` plotting extent.
#' @param degree_labels Use hemisphere degree labels automatically for
#'   longitude-latitude rasters, always, or never.
#' @param study_region Optional study-region boundary supplied as an `sf`,
#'   `sfc`, `Spatial`, `SpatVector`, or data frame with `x` and `y` columns.
#' @param region_colour,region_linewidth,region_linetype Appearance of the
#'   optional study-region boundary.
#'
#' @return A ggplot object.
#' @noRd
.plot_climniche_map <- function(x,
                                metric = c("niche_distance_change",
                                           "niche_boundary_exceedance",
                                           "climate_change_amount",
                                           "climate_reconfiguration",
                                           "change_alignment",
                                           "outside_niche_exceedance",
                                           "composition_change"),
                                occupied = NULL,
                                occupied_only = FALSE,
                                occupied_threshold = NULL,
                                title = NULL,
                                midpoint = 0,
                                limits = NULL,
                                breaks = NULL,
                                colours = NULL,
                                legend_title = FALSE,
                                legend_position = "right",
                                show_legend = TRUE,
                                symmetric = NULL,
                                extent = NULL,
                                degree_labels = c("auto", "none", "hemisphere"),
                                study_region = NULL,
                                region_colour = "black",
                                region_linewidth = 0.35,
                                region_linetype = 1) {
  .need_ggplot2()
  metric <- match.arg(metric)
  metric_key <- .metric_key(metric)

  if (inherits(x, "climniche_fit")) {
    if (is.null(x$rasters)) {
      stop("x does not contain raster outputs.", call. = FALSE)
    }
    r <- x$rasters[[metric_key]]
    if (is.null(r)) {
      r <- x$rasters[[metric]]
    }
  } else {
    r <- x
  }
  if (!.is_raster(r) && !.is_spatraster(r)) {
    stop("x must be a climniche object with raster outputs, a RasterLayer or a SpatRaster.",
         call. = FALSE)
  }
  if (occupied_only) {
    if (is.null(occupied)) {
      stop("occupied must be supplied when occupied_only = TRUE.",
           call. = FALSE)
    }
    r <- .mask_to_occupied(r, occupied,
                           occupied_threshold = occupied_threshold)
  }

  dat <- .spatial_df(r, value = "value")
  occ <- .occupied_df(occupied, occupied_threshold = occupied_threshold)
  lim <- .map_extent(dat, extent = extent)
  ttl <- .plot_title(title, .metric_label(metric))
  cell_size <- .spatial_res(r)
  axis_spec <- .map_axis_spec(r, degree_labels = degree_labels)
  region_df <- .study_region_df(study_region, r)
  if (isFALSE(legend_title)) {
    legend_title <- NULL
  } else if (is.null(legend_title)) {
    legend_title <- .metric_label(metric)
  }
  if (is.null(colours)) {
    colours <- .metric_map_colours(metric_key)
  }
  if (is.null(symmetric)) {
    symmetric <- metric_key %in% c("niche_distance_change", "change_alignment")
  }
  if (is.null(limits)) {
    if (identical(metric_key, "change_alignment")) {
      limits <- c(-1, 1)
    } else if (isTRUE(symmetric)) {
      radius <- max(abs(dat$value - midpoint), na.rm = TRUE)
      if (!is.finite(radius) || radius <= 0) radius <- 1
      limits <- midpoint + c(-radius, radius)
    } else if (metric_key %in% c("climate_change_amount",
                                 "climate_reconfiguration",
                                 "niche_boundary_exceedance")) {
      upper <- max(dat$value, na.rm = TRUE)
      if (!is.finite(upper) || upper <= 0) upper <- 1
      limits <- c(0, upper)
    }
  }
  scale_breaks <- if (is.null(breaks)) ggplot2::waiver() else breaks

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile(width = cell_size[1], height = cell_size[2]) +
    ggplot2::coord_equal(xlim = lim$xlim, ylim = lim$ylim, expand = FALSE) +
    ggplot2::scale_x_continuous(labels = axis_spec$x) +
    ggplot2::scale_y_continuous(labels = axis_spec$y) +
    ggplot2::labs(x = axis_spec$xlab, y = axis_spec$ylab, fill = legend_title,
                  title = ttl) +
    .map_theme() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) legend_position else "none"
    )

  if (metric_key %in% c("niche_distance_change", "change_alignment")) {
    p <- p + ggplot2::scale_fill_gradient2(
      low = colours[1], mid = colours[2], high = colours[3],
      midpoint = midpoint, limits = limits, breaks = scale_breaks,
      na.value = NA,
      oob = .squish_values
    )
  } else {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = colours, limits = limits, breaks = scale_breaks,
      na.value = NA,
      oob = .squish_values
    )
  }

  if (!occupied_only && !is.null(occ) && nrow(occ) > 0) {
    p <- p + ggplot2::geom_point(
      data = occ,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      size = 0.05,
      colour = "grey10",
      alpha = 0.14
    )
  }
  if (!is.null(region_df) && nrow(region_df) > 0) {
    p <- p + ggplot2::geom_path(
      data = region_df,
      ggplot2::aes(x = x, y = y, group = group),
      inherit.aes = FALSE,
      colour = region_colour,
      linewidth = region_linewidth,
      linetype = region_linetype,
      lineend = "round"
    )
  }
  p
}

#' Plot the climniche exposure plane
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells.
#' @param max_points Maximum number of points to draw.
#' @param seed Random seed used when subsampling.
#' @param title Optional plot title.
#' @param colour_by Colour cells by Niche Boundary Exceedance.
#'
#' @return A ggplot object.
#' @noRd
.plot_climniche_exposure <- function(x, scope = c("current", "all"),
                                     max_points = 6000, seed = 1,
                                     title = NULL,
                                     colour_by = "niche_boundary_exceedance") {
  .need_ggplot2()
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  dat <- .climniche_plot_data(x, scope = scope, max_points = max_points,
                              seed = seed)
  if (!identical(colour_by, "niche_boundary_exceedance")) {
    stop("colour_by must be 'niche_boundary_exceedance'.", call. = FALSE)
  }
  dat$x_value <- dat$climate_change_amount
  dat$y_value <- dat$niche_distance_change
  ttl <- .plot_title(title, "Niche relative climate exposure")

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x_value, y = y_value)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45",
                        linewidth = 0.3, linetype = 2) +
    ggplot2::geom_point(
      ggplot2::aes(colour = niche_boundary_exceedance),
      alpha = 0.42,
      size = 0.65
    ) +
    ggplot2::labs(
      x = "Climatic Displacement",
      y = "Niche Distance Shift",
      title = ttl,
      subtitle = "Positive Niche Distance Shift indicates movement away from the realised niche centre"
    ) +
    .climniche_theme()
  p <- p + ggplot2::scale_colour_gradientn(
    colours = .metric_map_colours("niche_boundary_exceedance"),
    name = "Niche Boundary\nExceedance"
  )
  p
}

#' Plot a climniche reported quantity distribution
#'
#' @param x A fitted `climniche_fit` object.
#' @param metric One of `"climate_change_amount"`, `"niche_distance_change"`,
#'   `"niche_boundary_exceedance"` or `"climate_reconfiguration"`.
#' @param scope `"current"` for a reference-weighted distribution or `"all"`
#'   for an unweighted distribution across evaluated cells.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#' @noRd
.plot_climniche_distribution <- function(x,
                                         metric = c("niche_distance_change",
                                                    "climate_change_amount",
                                                    "niche_boundary_exceedance",
                                                    "climate_reconfiguration",
                                                    "outside_niche_exceedance",
                                                    "composition_change"),
                                         scope = c("current", "all"),
                                         title = NULL) {
  .need_ggplot2()
  metric <- match.arg(metric)
  metric_key <- .metric_key(metric)
  scope <- match.arg(scope)
  tab <- climniche_table(x, scope = scope)
  tab$value <- tab[[metric_key]]
  tab$plot_weight <- if (scope == "current") {
    tab$occupied_weight
  } else {
    rep(1, nrow(tab))
  }
  ttl <- .plot_title(title, .metric_label(metric))
  y_label <- if (scope == "current") "Weighted count" else "Number of cells"

  ggplot2::ggplot(tab, ggplot2::aes(x = value, weight = plot_weight)) +
    ggplot2::geom_histogram(bins = 36, fill = "#5b91bd",
                            colour = "white", linewidth = 0.15) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey35",
                        linetype = 2, linewidth = 0.3) +
    ggplot2::labs(x = .metric_label(metric), y = y_label,
                  title = ttl) +
    .climniche_theme() +
    ggplot2::theme(legend.position = "none")
}

#' Plot a climniche report figure
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells.
#'
#' @return A patchwork object when `patchwork` is installed, otherwise a named
#'   list of ggplot objects.
#' @noRd
.plot_climniche_report <- function(x, scope = c("current", "all")) {
  scope <- match.arg(scope)
  plot_climniche_summary_figure(x, scope = scope)
}

#' Plot the four reported climniche quantities as maps
#'
#' @param x A fitted `climniche_fit` object with raster outputs.
#' @param metrics Character vector of reported quantity field names.
#' @param ncol Number of map columns when patchwork is available.
#' @param legend_title Shared legend title. The default suppresses repeated
#'   metric names because each panel is titled.
#' @param ... Map options passed to `plot_climniche_map()`. Commonly used
#'   arguments are `occupied`, `occupied_only`, `occupied_threshold`, `extent`,
#'   `degree_labels`, `study_region`, `legend_position`, `limits`, `breaks`,
#'   `colours`, and `show_legend`.
#'
#' @return A patchwork object when patchwork is installed, otherwise a named
#'   list of ggplot objects.
#' @export
plot_climniche_maps <- function(
    x,
    metrics = c("climate_change_amount", "niche_distance_change",
                "climate_reconfiguration", "niche_boundary_exceedance"),
    ncol = 2L,
    legend_title = FALSE,
    ...) {
  titles <- c(
    climate_change_amount = "Climatic Displacement",
    niche_distance_change = "Niche Distance Shift",
    climate_reconfiguration = "Climatic Reconfiguration",
    niche_boundary_exceedance = "Niche Boundary Exceedance"
  )
  metrics <- vapply(metrics, .metric_key, character(1))
  if (!length(metrics) || any(!metrics %in% names(titles))) {
    stop("metrics must contain reported metric field names or their legacy aliases.",
         call. = FALSE)
  }
  ncol <- .check_finite_scalar(ncol, "ncol")
  if (ncol != floor(ncol) || ncol < 1) {
    stop("ncol must be a positive integer.", call. = FALSE)
  }
  plots <- lapply(metrics, function(metric) {
    plot_climniche_map(
      x,
      metric = metric,
      title = unname(titles[metric]),
      legend_title = legend_title,
      ...
    )
  })
  names(plots) <- metrics
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(plots)
  }
  patchwork::wrap_plots(plots, ncol = as.integer(ncol)) +
    patchwork::plot_annotation(
      tag_levels = "a", tag_prefix = "(", tag_suffix = ")"
    )
}

#' Plot mean climatic variable contribution
#'
#' @param x A `climniche_fit` object.
#' @param occupied_only If TRUE, summarize occupied cells only.
#' @param variable_labels Optional named vector replacing variable labels.
#' @param title Optional plot title. Use `FALSE` to suppress it.
#'
#' @return A ggplot object.
#' @export
plot_climniche_variable_contribution <- function(x, occupied_only = TRUE,
                                                 variable_labels = NULL,
                                                 title = NULL) {
  .need_ggplot2()
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a climniche object.", call. = FALSE)
  }
  idx <- if (occupied_only) x$occupied else seq_len(nrow(x$variable_contribution))
  weights <- if (occupied_only) {
    .fit_reference_weights(x)[idx]
  } else {
    rep(1, length(idx))
  }
  vals <- .weighted_col_means(x$variable_contribution[idx, , drop = FALSE],
                              weights)
  variables <- names(vals)
  if (!is.null(variable_labels)) {
    replace <- match(variables, names(variable_labels))
    hit <- !is.na(replace)
    variables[hit] <- unname(variable_labels[replace[hit]])
  }
  dat <- data.frame(
    variable = factor(variables, levels = variables),
    contribution = as.numeric(vals)
  )

  ggplot2::ggplot(dat, ggplot2::aes(x = variable, y = contribution,
                                    fill = contribution > 0)) +
    ggplot2::geom_col(width = 0.65, colour = "grey25", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#c65d57",
                                          "FALSE" = "#4c78a8")) +
    ggplot2::labs(
      x = NULL,
      y = "Mean contribution to niche potential change",
      title = .plot_title(title, "Contributions to niche potential change")
    ) +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(colour = "black"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                          colour = "black", size = 7),
      axis.text.y = ggplot2::element_text(colour = "black", size = 7),
      axis.title.y = ggplot2::element_text(colour = "black", size = 8),
      axis.line = ggplot2::element_line(linewidth = 0.2, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.2, colour = "black")
    )
}
