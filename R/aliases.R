#' Fit niche climate exposure
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
#' @param boundary Quantile defining the empirical realised niche boundary.
#' @param scale Logical. If TRUE, standardise current and future values.
#' @param global_mean Optional means used for standardisation.
#' @param global_sd Optional standard deviations used for standardisation.
#' @param tolerance Optional Niche Distance Shift tolerance.
#' @param tolerance_quantile Quantile of absolute Niche Distance Shift used
#'   when `tolerance = NULL`.
#' @param stable_climate_change Optional threshold for limited climate niche
#'   change.
#' @param stable_quantile Quantile of Climatic Displacement used when
#'   `stable_climate_change = NULL`.
#' @param stable_reconfiguration Optional threshold for low Climatic
#'   Reconfiguration.
#' @param stable_reconfiguration_quantile Quantile of Climatic Reconfiguration
#'   used when `stable_reconfiguration = NULL`.
#' @param boundary_exceedance_tolerance Tolerance for deciding whether future
#'   climate exceeds the empirical niche boundary.
#' @param conflict_ratio Minimum minority sign contribution share used to mark
#'   mixed variable responses. Set to NULL to disable this flag.
#'
#' @return An object of class `climniche_fit`.
#'
#' @details
#' The fitted object stores four primary metrics as snake_case fields:
#' Climatic Displacement (`climate_change_amount`), Niche Distance Shift
#' (`niche_distance_change`), Climatic Reconfiguration
#' (`climate_reconfiguration`) and Niche Boundary Exceedance
#' (`niche_boundary_exceedance`). Let current and future climatic conditions at
#' cell `i` be `c_i` and `f_i`, let `mu` be the centre of the current realised
#' climatic niche, and let `d_A(x, y)` be the sensitivity weighted distance
#' under weighting matrix `A`. Climatic Displacement is `d_A(f_i, c_i)`. Niche
#' Distance Shift is `d_A(f_i, mu) - d_A(c_i, mu)`. Climatic Reconfiguration is
#' `sqrt(max(0, D_i^2 - R_i^2))`, where `D_i` is Climatic Displacement and
#' `R_i` is Niche Distance Shift. Niche Boundary Exceedance is
#' `max(0, d_A(f_i, mu) - B_q)`, where `B_q` is the `q`-th weighted quantile of
#' current reference cell distances from the realised niche centre.
#' @export
fit_climniche <- function(current, future, occupied = NULL,
                          occupied_threshold = NULL, cnfa = NULL,
                          center = NULL, sensitivity = NULL, A = NULL,
                          metric = c("diag", "factor"),
                          boundary = 0.95, scale = TRUE,
                          global_mean = NULL, global_sd = NULL,
                          tolerance = NULL,
                          tolerance_quantile = 0.10,
                          stable_climate_change = NULL,
                          stable_quantile = 0.25,
                          stable_reconfiguration = NULL,
                          stable_reconfiguration_quantile = 0.25,
                          boundary_exceedance_tolerance = 0,
                          conflict_ratio = 0.25) {
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
    stable_climate_change = stable_climate_change,
    stable_quantile = stable_quantile,
    stable_reconfiguration = stable_reconfiguration,
    stable_reconfiguration_quantile = stable_reconfiguration_quantile,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance,
    conflict_ratio = conflict_ratio
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
#' @param ... Additional arguments passed to `fit_climniche()`.
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
#' @param ... Additional arguments passed to `fit_climniche()`.
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
#' @param x A fitted climniche object with raster outputs, or a RasterLayer.
#' @param metric Metric to plot.
#' @param occupied Optional current reference RasterLayer to overlay.
#' @param occupied_only If TRUE, mask the plotted raster to current occurrence
#'   cells.
#' @param occupied_threshold Threshold used when `occupied` contains binary or
#'   continuous values. Values above the threshold keep their original value
#'   when used as an overlay or mask.
#' @param title Optional plot title. Use `FALSE` to suppress it.
#' @param midpoint Midpoint for the Niche Distance Shift colour scale.
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
                               midpoint = 0) {
  metric <- match.arg(metric)
  .plot_climniche_map(x = x, metric = metric, occupied = occupied,
                      occupied_only = occupied_only,
                      occupied_threshold = occupied_threshold,
                      title = title,
                      midpoint = midpoint)
}

#' Plot climniche classes
#'
#' @param x A fitted climniche object with raster outputs.
#' @param occupied Optional current reference RasterLayer to overlay.
#' @param occupied_only If TRUE, mask the plotted classes to current occurrence
#'   cells.
#' @param occupied_threshold Threshold used when `occupied` contains binary or
#'   continuous values. Values above the threshold keep their original value
#'   when used as an overlay or mask.
#' @param title Optional plot title. Use `FALSE` to suppress it.
#'
#' @return A ggplot object.
#' @export
plot_climniche_classes <- function(x, occupied = NULL, occupied_only = FALSE,
                                   occupied_threshold = NULL,
                                   title = NULL) {
  .plot_climniche_classes(x = x, occupied = occupied,
                          occupied_only = occupied_only,
                          occupied_threshold = occupied_threshold,
                          title = title)
}

#' Plot the climniche exposure plane
#'
#' @param x A fitted climniche object.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells.
#' @param max_points Maximum number of points to draw.
#' @param seed Random seed used when subsampling.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#' @export
plot_climniche_exposure <- function(x, scope = c("current", "all"),
                                    max_points = 6000, seed = 1,
                                    title = NULL) {
  .plot_climniche_exposure(x = x, scope = scope, max_points = max_points,
                           seed = seed, title = title)
}

#' Plot climniche class proportions
#'
#' @param x A fitted climniche object.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#' @export
plot_climniche_class_summary <- function(x, scope = c("current", "all"),
                                         title = NULL) {
  .plot_climniche_class_summary(x = x, scope = scope, title = title)
}

#' Plot a climniche metric distribution
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
