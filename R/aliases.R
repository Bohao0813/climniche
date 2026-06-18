#' Fit niche relative climate exposure
#'
#' @param current Numeric matrix or data frame of current environmental values.
#' @param future Numeric matrix or data frame of future environmental values.
#' @param occupied NULL, logical vector, row indices, or a numeric vector with
#'   one value per row identifying current occurrence, range, or continuous SDM
#'   suitability weights.
#' @param occupied_threshold Threshold used when `occupied` contains binary or
#'   continuous values. Values at or below the threshold receive zero reference
#'   weight; values above it keep their original value.
#' @param cnfa Optional CENFA model object.
#' @param center Optional realised niche centre in standardised climate space.
#' @param sensitivity Optional environmental sensitivity weights.
#' @param A Optional niche metric matrix.
#' @param metric Metric construction when `A` is not supplied.
#' @param boundary User-set quantile defining the empirical realised niche
#'   boundary. The default is `0.95`.
#' @param scale Logical. If TRUE, standardise current and future values.
#' @param global_mean Optional means used for standardisation.
#' @param global_sd Optional standard deviations used for standardisation.
#' @param tolerance Optional user-set Niche Distance Shift tolerance. If `NULL`,
#'   the value is calculated from `tolerance_quantile`.
#' @param tolerance_quantile Quantile of absolute Niche Distance Shift used
#'   when `tolerance = NULL`.
#' @param boundary_exceedance_tolerance User-set tolerance for deciding whether
#'   future climate exceeds the empirical niche boundary.
#'
#' @return An object of class `climniche_fit`.
#'
#' @details
#' The fitted object stores four related quantities as snake_case fields:
#' Climatic Displacement (`climate_change_amount`), Niche Distance Shift
#' (`niche_distance_change`), Climatic Reconfiguration
#' (`climate_reconfiguration`) and Niche Boundary Exceedance
#' (`niche_boundary_exceedance`).
#'
#' Let current and future climatic conditions at cell \eqn{i} be \eqn{c_i} and
#' \eqn{f_i}. Let \eqn{\mu} be the centre of the current realised climatic niche,
#' and let \eqn{d_A(x, y)} be the sensitivity weighted distance under weighting
#' matrix \eqn{A}. The niche-relative decomposition reports
#'
#' \deqn{D_i = d_A(f_i, c_i)}
#'
#' \deqn{R_i = d_A(f_i, \mu) - d_A(c_i, \mu)}
#'
#' \deqn{C_i = \sqrt{\max(0, D_i^2 - R_i^2)}}
#'
#' \deqn{E_i = \max(0, d_A(f_i, \mu) - B_q)}
#'
#' where \eqn{B_q} is the \eqn{q}-th weighted quantile of current reference cell
#' distances from the realised niche centre. Positive Niche Boundary Exceedance
#' is therefore an excess distance beyond this empirical radial boundary.
#' Climatic Reconfiguration is derived from Climatic Displacement and Niche
#' Distance Shift; it is a non radial displacement component rather than an
#' independently estimated process.
#'
#' Descriptor thresholds are user-settable. If `tolerance = NULL`,
#' `climniche` calculates the effective Niche Distance Shift tolerance from
#' `tolerance_quantile`. The fitted object stores the effective values in
#' `descriptor_settings`.
#'
#' @section User-settable thresholds:
#' - `boundary`: quantile used to define the empirical realised niche boundary.
#' - `tolerance`: direct Niche Distance Shift tolerance; otherwise
#'   `tolerance_quantile`.
#' - `boundary_exceedance_tolerance`: direct tolerance for Niche Boundary
#'   Exceedance.
#' @export
fit_climniche <- function(current, future, occupied = NULL,
                          occupied_threshold = NULL, cnfa = NULL,
                          center = NULL, sensitivity = NULL, A = NULL,
                          metric = c("diag", "factor"),
                          boundary = 0.95, scale = TRUE,
                          global_mean = NULL, global_sd = NULL,
                          tolerance = NULL,
                          tolerance_quantile = 0.10,
                          boundary_exceedance_tolerance = 0) {
  .fit_climniche_matrix(
    current = current,
    future = future,
    occupied = occupied,
    occupied_threshold = occupied_threshold,
    cnfa = cnfa,
    center = center,
    sensitivity = sensitivity,
    A = A,
    metric = metric,
    boundary = boundary,
    scale = scale,
    global_mean = global_mean,
    global_sd = global_sd,
    tolerance = tolerance,
    tolerance_quantile = tolerance_quantile,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance
  )
}

#' Fit climniche to raster data
#'
#' @param current Raster* object of current environmental layers.
#' @param future Raster* object of future environmental layers.
#' @param occupied Optional RasterLayer with binary or continuous occurrence,
#'   range, or SDM suitability values.
#' @param occupied_threshold Values at or below this threshold receive zero
#'   reference weight. Values above it keep their original value.
#' @param domain Optional RasterLayer limiting cells where exposure is analysed.
#' @param domain_threshold Values greater than this threshold define the domain.
#' @param ... Additional arguments passed to `fit_climniche()`, including
#'   `boundary`, `tolerance`, `tolerance_quantile`, and
#'   `boundary_exceedance_tolerance`.
#'
#' @return An object of class `climniche_fit` with raster outputs.
#' @export
fit_climniche_raster <- function(current, future, occupied = NULL,
                                 occupied_threshold = NULL, domain = NULL,
                                 domain_threshold = 0, ...) {
  .fit_climniche_raster(current = current, future = future,
                        occupied = occupied,
                        occupied_threshold = occupied_threshold,
                        domain = domain,
                        domain_threshold = domain_threshold, ...)
}

#' Fit climniche to terra raster data
#'
#' @param current terra SpatRaster of current environmental layers.
#' @param future terra SpatRaster of future environmental layers.
#' @param occupied Optional one layer SpatRaster with binary or continuous
#'   occurrence, range, or SDM suitability values.
#' @param occupied_threshold Values at or below this threshold receive zero
#'   reference weight. Values above it keep their original value.
#' @param domain Optional one layer SpatRaster limiting cells where exposure is
#'   analysed.
#' @param domain_threshold Values greater than this threshold define the domain.
#' @param ... Additional arguments passed to `fit_climniche()`, including
#'   `boundary`, `tolerance`, `tolerance_quantile`, and
#'   `boundary_exceedance_tolerance`.
#'
#' @return An object of class `climniche_fit` with raster outputs.
#' @export
fit_climniche_terra <- function(current, future, occupied = NULL,
                                occupied_threshold = NULL, domain = NULL,
                                domain_threshold = 0, ...) {
  .fit_climniche_terra(current = current, future = future,
                       occupied = occupied,
                       occupied_threshold = occupied_threshold,
                       domain = domain,
                       domain_threshold = domain_threshold, ...)
}

#' Plot a climniche map
#'
#' @param x A fitted climniche object with raster outputs, a RasterLayer, or a
#'   terra SpatRaster.
#' @param metric Metric to plot.
#' @param occupied Optional current reference RasterLayer or terra SpatRaster
#'   to overlay.
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
#' @param legend_title Optional colour legend title.
#' @param legend_position Position passed to the ggplot theme.
#' @param show_legend If FALSE, suppress the colour legend.
#' @param symmetric If TRUE, use limits symmetric around `midpoint`.
#' @param extent Optional `c(xmin, xmax, ymin, ymax)` plotting extent.
#' @param degree_labels Use hemisphere degree labels automatically for
#'   longitude-latitude rasters, always, or never.
#' @param study_region Optional study-region boundary supplied as an `sf`,
#'   `sfc`, `Spatial`, `SpatVector`, or data frame with `x` and `y` columns.
#' @param region_colour,region_linewidth,region_linetype Appearance of the
#'   optional study-region boundary.
#'
#' @return A ggplot object.
#' @export
plot_climniche_map <- function(x,
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
                               legend_title = NULL,
                               legend_position = "right",
                               show_legend = TRUE,
                               symmetric = NULL,
                               extent = NULL,
                               degree_labels = c("auto", "none", "hemisphere"),
                               study_region = NULL,
                               region_colour = "black",
                               region_linewidth = 0.35,
                               region_linetype = 1) {
  metric <- match.arg(metric)
  .plot_climniche_map(x = x, metric = metric, occupied = occupied,
                      occupied_only = occupied_only,
                      occupied_threshold = occupied_threshold,
                      title = title,
                      midpoint = midpoint,
                      limits = limits,
                      breaks = breaks,
                      colours = colours,
                      legend_title = legend_title,
                      legend_position = legend_position,
                      show_legend = show_legend,
                      symmetric = symmetric,
                      extent = extent,
                      degree_labels = degree_labels,
                      study_region = study_region,
                      region_colour = region_colour,
                      region_linewidth = region_linewidth,
                      region_linetype = region_linetype)
}

#' Plot the climniche exposure plane
#'
#' @param x A fitted climniche object.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells.
#' @param max_points Maximum number of points to draw.
#' @param seed Random seed used when subsampling.
#' @param title Optional plot title.
#' @param colour_by Colour cells by Niche Boundary Exceedance.
#'
#' @return A ggplot object.
#' @export
plot_climniche_exposure <- function(x, scope = c("current", "all"),
                                    max_points = 6000, seed = 1,
                                    title = NULL,
                                    colour_by = "niche_boundary_exceedance") {
  .plot_climniche_exposure(x = x, scope = scope, max_points = max_points,
                           seed = seed, title = title,
                           colour_by = colour_by)
}

#' Plot a climniche reported quantity distribution
#'
#' @param x A fitted climniche object.
#' @param metric Metric to plot.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#' @export
plot_climniche_distribution <- function(x,
                                        metric = c("niche_distance_change",
                                                   "climate_change_amount",
                                                   "niche_boundary_exceedance",
                                                   "climate_reconfiguration",
                                                   "outside_niche_exceedance",
                                                   "composition_change"),
                                        scope = c("current", "all"),
                                        title = NULL) {
  metric <- match.arg(metric)
  .plot_climniche_distribution(x = x, metric = metric, scope = scope,
                               title = title)
}

#' Plot a climniche report figure
#'
#' @param x A fitted climniche object.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells.
#'
#' @return A patchwork object when `patchwork` is installed, otherwise a named
#'   list of ggplot objects.
#' @export
plot_climniche_report <- function(x, scope = c("current", "all")) {
  .plot_climniche_report(x = x, scope = scope)
}
