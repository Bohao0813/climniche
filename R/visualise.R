.metric_label <- function(metric) {
  switch(metric,
         climate_change_amount = "Climate change amount",
         niche_distance_change = "Niche distance change",
         composition_change = "Composition change",
         change_alignment = "Change alignment",
         outside_niche_exceedance = "Niche boundary exceedance",
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

.raster_df <- function(x, value = "value") {
  .need_raster()
  pts <- raster::rasterToPoints(x)
  pts <- as.data.frame(pts)
  names(pts) <- c("x", "y", value)
  pts
}

.occupied_df <- function(occupied) {
  if (is.null(occupied)) {
    return(NULL)
  }
  pts <- .raster_df(occupied, value = "occupied")
  pts[!is.na(pts$occupied) & pts$occupied > 0, , drop = FALSE]
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

.class_colours <- function() {
  c(
    "little climate niche change" = "#f0f0f0",
    "closer to current niche" = "#4c78a8",
    "farther from current niche" = "#c65d57",
    "outside current niche boundary" = "#d99b45",
    "changed composition, similar distance" = "#5f9e8f"
  )
}

.climniche_theme <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 0.2),
      plot.subtitle = ggplot2::element_text(colour = "grey35", size = base_size - 0.8),
      axis.text = ggplot2::element_text(colour = "grey30", size = base_size - 1.4),
      axis.title = ggplot2::element_text(colour = "grey15", size = base_size - 0.5),
      axis.line = ggplot2::element_line(linewidth = 0.25),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = base_size - 1.5),
      legend.title = ggplot2::element_text(size = base_size - 1.2),
      legend.key.height = grid::unit(4.6, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.4)
    )
}

.map_theme <- function() {
  ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 8.7),
      axis.text = ggplot2::element_text(colour = "grey35", size = 6.4),
      axis.ticks = ggplot2::element_line(colour = "grey70", linewidth = 0.2),
      axis.line = ggplot2::element_line(linewidth = 0.2),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 6.8),
      legend.title = ggplot2::element_text(size = 7.2),
      legend.key.width = grid::unit(3.4, "mm"),
      legend.key.height = grid::unit(8, "mm")
    )
}

.climniche_plot_data <- function(x, scope = c("current", "all"),
                                 max_points = 5000, seed = 1) {
  scope <- match.arg(scope)
  dat <- climniche_table(x, scope = scope)
  if (nrow(dat) > max_points) {
    set.seed(seed)
    dat <- dat[sample(seq_len(nrow(dat)), max_points), , drop = FALSE]
  }
  dat$class <- factor(dat$class, levels = names(.class_colours()))
  dat
}

#' Plot a climniche raster metric
#'
#' @param x A `climniche_fit` object with raster outputs, or a RasterLayer.
#' @param metric Metric to plot when `x` is a `climniche_fit` object.
#' @param occupied Optional current occurrence/range RasterLayer to overlay.
#' @param occupied_only If TRUE, mask the plotted raster to current occurrence
#'   cells.
#' @param title Optional plot title. Use `FALSE` to suppress it.
#' @param midpoint Midpoint for the niche distance change colour scale.
#'
#' @return A ggplot object.
.plot_climniche_map <- function(x,
                                metric = c("niche_distance_change",
                                           "outside_niche_exceedance",
                                           "climate_change_amount",
                                           "composition_change",
                                           "change_alignment"),
                                occupied = NULL,
                                occupied_only = FALSE,
                                title = NULL,
                                midpoint = 0) {
  .need_ggplot2()
  .need_raster()
  metric <- match.arg(metric)

  if (inherits(x, "climniche_fit")) {
    if (is.null(x$rasters)) {
      stop("x does not contain raster outputs.", call. = FALSE)
    }
    r <- x$rasters[[metric]]
  } else {
    r <- x
  }
  if (!methods::is(r, "Raster")) {
    stop("x must be a climniche object with rasters, or a RasterLayer.",
         call. = FALSE)
  }
  if (occupied_only) {
    if (is.null(occupied)) {
      stop("occupied must be supplied when occupied_only = TRUE.",
           call. = FALSE)
    }
    r <- raster::mask(r, occupied)
  }

  dat <- .raster_df(r, value = "value")
  occ <- .occupied_df(occupied)
  lim <- .crop_df(dat)
  ttl <- .plot_title(title, .metric_label(metric))
  cell_size <- raster::res(r)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile(width = cell_size[1], height = cell_size[2]) +
    ggplot2::coord_equal(xlim = lim$xlim, ylim = lim$ylim, expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL, title = ttl) +
    .map_theme()

  if (identical(metric, "niche_distance_change")) {
    p <- p + ggplot2::scale_fill_gradient2(
      low = "#4c78a8", mid = "white", high = "#c65d57",
      midpoint = midpoint, na.value = NA
    )
  } else if (identical(metric, "change_alignment")) {
    p <- p + ggplot2::scale_fill_gradient2(
      low = "#4c78a8", mid = "white", high = "#c65d57",
      midpoint = 0, limits = c(-1, 1), na.value = NA
    )
  } else if (identical(metric, "outside_niche_exceedance")) {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = c("white", "#f2d7a0", "#d99b45", "#7a3b20"),
      na.value = NA
    )
  } else if (identical(metric, "composition_change")) {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = c("#f7fcf5", "#b7e0c2", "#5f9e8f", "#1b5e50"),
      na.value = NA
    )
  } else {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = c("#f7fbff", "#aac7df", "#5b91bd", "#1f4f78"),
      na.value = NA
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
  p
}

#' Plot the climniche exposure plane
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current occurrence/range cells or `"all"` for
#'   all evaluated cells.
#' @param max_points Maximum number of points to draw.
#' @param seed Random seed used when subsampling.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
.plot_climniche_exposure <- function(x, scope = c("current", "all"),
                                     max_points = 6000, seed = 1,
                                     title = NULL) {
  .need_ggplot2()
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  dat <- .climniche_plot_data(x, scope = scope, max_points = max_points,
                              seed = seed)
  dat$x_value <- dat$climate_change_amount
  dat$y_value <- dat$niche_distance_change
  ttl <- .plot_title(title, "Niche climate exposure")

  ggplot2::ggplot(
    dat,
    ggplot2::aes(x = x_value, y = y_value, colour = class)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45",
                        linewidth = 0.3, linetype = 2) +
    ggplot2::geom_point(alpha = 0.42, size = 0.65) +
    ggplot2::scale_colour_manual(values = .class_colours(), drop = FALSE) +
    ggplot2::labs(
      x = "Climate change amount",
      y = "Niche distance change",
      colour = NULL,
      title = ttl,
      subtitle = "Positive values indicate divergence from the current realised climatic niche"
    ) +
    .climniche_theme()
}

#' Plot climniche class proportions
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current occurrence/range cells or `"all"` for
#'   all evaluated cells.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
.plot_climniche_class_summary <- function(x, scope = c("current", "all"),
                                          title = NULL) {
  .need_ggplot2()
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  tab <- climniche_table(x, scope = scope)
  dat <- as.data.frame(prop.table(table(tab$class)), stringsAsFactors = FALSE)
  names(dat) <- c("class", "proportion")
  dat$class <- factor(dat$class, levels = names(.class_colours()))
  dat <- dat[order(dat$proportion), , drop = FALSE]
  dat$prop <- dat$proportion
  ttl <- .plot_title(title, "climniche class proportions")

  ggplot2::ggplot(dat, ggplot2::aes(x = class, y = prop, fill = class)) +
    ggplot2::geom_col(width = 0.68, colour = "grey25", linewidth = 0.18) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = function(z) paste0(round(100 * z), "%")) +
    ggplot2::scale_fill_manual(values = .class_colours(), drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Proportion of analysed cells", title = ttl) +
    .climniche_theme() +
    ggplot2::theme(legend.position = "none")
}

#' Plot a climniche metric distribution
#'
#' @param x A fitted `climniche_fit` object.
#' @param metric One of `"climate_change_amount"`, `"niche_distance_change"`,
#'   `"outside_niche_exceedance"` or `"composition_change"`.
#' @param scope `"current"` for current occurrence/range cells or `"all"` for
#'   all evaluated cells.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
.plot_climniche_distribution <- function(x,
                                         metric = c("niche_distance_change",
                                                    "climate_change_amount",
                                                    "outside_niche_exceedance",
                                                    "composition_change"),
                                         scope = c("current", "all"),
                                         title = NULL) {
  .need_ggplot2()
  metric <- match.arg(metric)
  scope <- match.arg(scope)
  tab <- climniche_table(x, scope = scope)
  tab$value <- tab[[metric]]
  ttl <- .plot_title(title, .metric_label(metric))

  ggplot2::ggplot(tab, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 36, fill = "#5b91bd",
                            colour = "white", linewidth = 0.15) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey35",
                        linetype = 2, linewidth = 0.3) +
    ggplot2::labs(x = .metric_label(metric), y = "Number of cells",
                  title = ttl) +
    .climniche_theme() +
    ggplot2::theme(legend.position = "none")
}

#' Plot a compact climniche report figure
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current occurrence/range cells or `"all"` for
#'   all evaluated cells.
#'
#' @return A patchwork object when `patchwork` is installed, otherwise a named
#'   list of ggplot objects.
.plot_climniche_report <- function(x, scope = c("current", "all")) {
  scope <- match.arg(scope)
  plot_climniche_showcase(x, scope = scope)
}

#' Plot projected climate change classes
#'
#' @param x A `climniche_fit` object with raster outputs.
#' @param occupied Optional occupied-cell RasterLayer to overlay.
#' @param occupied_only If TRUE, mask the plotted classes to occupied cells.
#' @param title Optional plot title. Use `FALSE` to suppress it.
#'
#' @return A ggplot object.
.plot_climniche_classes <- function(x, occupied = NULL, occupied_only = FALSE,
                                    title = NULL) {
  .need_ggplot2()
  .need_raster()
  if (!inherits(x, "climniche_fit") || is.null(x$rasters$classification)) {
    stop("x must be a climniche object with a classification raster.",
         call. = FALSE)
  }
  r <- x$rasters$classification
  if (occupied_only) {
    if (is.null(occupied)) {
      stop("occupied must be supplied when occupied_only = TRUE.",
           call. = FALSE)
    }
    r <- raster::mask(r, occupied)
  }
  dat <- .raster_df(r, value = "class_id")
  lookup <- x$class_lookup
  dat$class <- lookup$class[match(dat$class_id, lookup$id)]
  dat$class <- factor(dat$class, levels = lookup$class)
  occ <- .occupied_df(occupied)
  lim <- .crop_df(dat)
  cell_size <- raster::res(r)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, fill = class)) +
    ggplot2::geom_tile(width = cell_size[1], height = cell_size[2]) +
    ggplot2::coord_equal(xlim = lim$xlim, ylim = lim$ylim, expand = FALSE) +
    ggplot2::scale_fill_manual(
      values = .class_colours(),
      na.value = NA
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL,
                  title = .plot_title(title, "Climate niche change classes")) +
    .map_theme() +
    ggplot2::theme(legend.key.height = grid::unit(4.6, "mm"))

  if (!occupied_only && !is.null(occ) && nrow(occ) > 0) {
    p <- p + ggplot2::geom_point(
      data = occ,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      size = 0.05,
      colour = "grey10",
      alpha = 0.12
    )
  }
  p
}

#' Plot mean variable contribution
#'
#' @param x A `climniche_fit` object.
#' @param occupied_only If TRUE, summarize occupied cells only.
#'
#' @return A ggplot object.
#' @export
plot_variable_contribution <- function(x, occupied_only = TRUE) {
  .need_ggplot2()
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a climniche object.", call. = FALSE)
  }
  idx <- if (occupied_only) x$occupied else seq_len(nrow(x$variable_contribution))
  vals <- colMeans(x$variable_contribution[idx, , drop = FALSE])
  dat <- data.frame(
    variable = factor(names(vals), levels = names(vals)),
    contribution = as.numeric(vals)
  )

  ggplot2::ggplot(dat, ggplot2::aes(x = variable, y = contribution,
                                    fill = contribution > 0)) +
    ggplot2::geom_col(width = 0.65, colour = "grey25", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#c65d57",
                                          "FALSE" = "#4c78a8")) +
    ggplot2::labs(x = NULL, y = "Mean contribution") +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.title.y = ggplot2::element_text(size = 8),
      axis.line = ggplot2::element_line(linewidth = 0.2),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )
}
