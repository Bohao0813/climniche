#' Fit niche relative climate exposure from matrices
#'
#' @param current Numeric matrix or data frame of current climate values. Rows
#'   are cells, sites, or samples. Columns are climate variables.
#' @param future Numeric matrix or data frame of future climate values with the
#'   same rows and variables as `current`. Complete row and variable names are
#'   matched before fitting.
#' @param occupied Reference information used to estimate the current realised
#'   niche. Use `NULL` to give every row weight 1, a logical vector to mark
#'   reference rows, a numeric vector of length `nrow(current)` for continuous
#'   reference weights, or positive integer row indices for 0/1 reference cells.
#'   Complete names are matched to the row names of `current`.
#' @param occupied_threshold Optional cutoff for numeric reference weights.
#'   Values at or below the cutoff are set to 0. Values above it keep their
#'   original continuous value.
#' @param cnfa Optional CENFA `cnfa` object. Its `mf` and `sf` components are
#'   used when `center`, `sensitivity`, or `A` are not supplied; `metric =
#'   "factor"` requires `co` and `eig`. The object and climatic inputs must use
#'   the same variables and standardisation.
#' @param center Optional realised niche centre on the scale used for distance
#'   calculations. With `scale = TRUE`, supply a centre in the standardized
#'   climatic space. If omitted, the centre is the weighted mean of current
#'   reference rows.
#' @param sensitivity Optional non-negative variable weights. Used to build a
#'   diagonal metric when `A` is not supplied.
#' @param A Optional square metric matrix defined for the fitted climatic
#'   space. When supplied, it overrides `sensitivity`, `cnfa`, and `metric` for
#'   distance calculations.
#' @param metric Method used to build `A` when `A` is missing. `"diag"` uses
#'   variable-level sensitivity weights. `"factor"` constructs a factor-weighted
#'   metric from the `co` and `eig` components of a compatible CENFA object.
#' @param boundary Quantile used to define the empirical boundary of the current
#'   realised niche. Must be between 0 and 1.
#' @param scale Logical. If `TRUE`, current and future values are centred and
#'   scaled with the current-layer mean and standard deviation before distances
#'   are calculated.
#' @param preprocess Logical. If `TRUE`, remove near-zero variance variables
#'   and highly correlated variables before metric fitting.
#' @param preprocess_correlation Maximum absolute pairwise correlation retained
#'   among current climate variables during preprocessing.
#' @param preprocess_min_sd Minimum current-climate standard deviation retained
#'   during preprocessing.
#' @param global_mean Optional means used for centering when `scale = TRUE`. If
#'   omitted, column means of `current` are used.
#' @param global_sd Optional standard deviations used for scaling. If
#'   `scale = TRUE` and this argument is omitted, column standard deviations of
#'   `current` are used.
#' @param tolerance Optional tolerance around zero for Niche Distance Shift. If
#'   `NULL`, the fitted object uses `tolerance_quantile`.
#' @param tolerance_quantile Quantile of absolute Niche Distance Shift used to
#'   set `tolerance` when `tolerance = NULL`.
#' @param boundary_exceedance_tolerance Tolerance used to label whether Niche
#'   Boundary Exceedance is greater than zero for descriptor summaries.
#'
#' @return A `climniche_fit` object containing the four metrics, reference
#'   weights, fitted niche centre, metric matrix and effective settings.
#'
#' @details
#' Let current and future climatic conditions at cell \eqn{i} be \eqn{c_i} and
#' \eqn{f_i}. Let \eqn{\mu} be the centre of the current realised climatic niche,
#' and let \eqn{d_A(x, y)} be the climatic distance under weighting matrix
#' \eqn{A}. The four metrics are:
#'
#' Climatic Displacement:
#' \deqn{D_i = d_A(f_i, c_i)}
#'
#' Niche Distance Shift:
#' \deqn{R_i = d_A(f_i, \mu) - d_A(c_i, \mu)}
#'
#' Climatic Reconfiguration:
#' \deqn{C_i = \sqrt{\max(0, D_i^2 - R_i^2)}}
#'
#' Niche Boundary Exceedance:
#' \deqn{E_i = \max(0, d_A(f_i, \mu) - B_q)}
#'
#' where \eqn{B_q} is the \eqn{q}-th weighted quantile of current reference cell
#' distances from the realised niche centre. It is the smallest observed
#' reference distance whose cumulative normalised reference weight reaches
#' \eqn{q}.
#'
#' In the transformed climatic space, let \eqn{r_{0i}} and \eqn{r_{1i}} be the
#' current and future distances from the niche centre, and let \eqn{\theta_i}
#' be the angle between their centred vectors. Then
#' \deqn{C_i^2 = 2 r_{0i} r_{1i} (1 - \cos(\theta_i)).}
#' Climatic Reconfiguration therefore combines angular change with current and
#' future niche distances. It is calculated from Climatic Displacement and
#' Niche Distance Shift rather than fitted independently.
#'
#' @section Metric fields:
#' The primary fitted fields are `climate_change_amount`,
#' `niche_distance_change`, `climate_reconfiguration`, and
#' `niche_boundary_exceedance`. The legacy names `composition_change` and
#' `outside_niche_exceedance` are retained as aliases for old code.
#'
#' @section User-settable thresholds:
#' `boundary` controls the empirical niche boundary. `tolerance` controls the
#' zero band for Niche Distance Shift. `boundary_exceedance_tolerance` controls
#' the boundary descriptor. The fitted values are stored in
#' `descriptor_settings`.
#'
#' @section Boundary scale:
#' `boundary_distance` and `boundary_radius` store the fitted boundary in
#' distance units. `boundary_potential` stores its squared value and is the
#' quantity accepted by [boundary_exceedance()]. `boundary_value` is retained
#' as a legacy alias of `boundary_distance`.
#'
#' @section Scaling and preprocessing:
#' `preprocess` selects retained variables. `scale` then converts those
#' variables to z scores using current-climate means and standard deviations.
#' Both are enabled by default.
#'
#' @section Choosing a fit function:
#' Use `fit_climniche()` when current and future climate values have already
#' been extracted to matrices or data frames. Use `fit_climniche_raster()` for
#' objects from the `raster` package. Use `fit_climniche_terra()` for
#' `terra::SpatRaster` objects. The three functions calculate the same metrics;
#' the raster and terra methods add spatial masking and return map layers in
#' `x$rasters`. Use [fit_climniche_series()] for ordered future periods or
#' climate model ensembles.
#'
#' @examples
#' sim <- simulate_climniche(n = 250, p = 6, seed = 7)
#' fit <- fit_climniche(
#'   current = sim$current,
#'   future = sim$future_away,
#'   occupied = sim$occupied,
#'   sensitivity = sim$sensitivity
#' )
#' climniche_summary(fit)
#' @export
fit_climniche <- function(current, future, occupied = NULL,
                          occupied_threshold = NULL, cnfa = NULL,
                          center = NULL, sensitivity = NULL, A = NULL,
                          metric = c("diag", "factor"),
                          boundary = 0.95, scale = TRUE,
                          preprocess = TRUE,
                          preprocess_correlation = 0.95,
                          preprocess_min_sd = 1e-08,
                          global_mean = NULL, global_sd = NULL,
                          tolerance = NULL,
                          tolerance_quantile = 0.10,
                          boundary_exceedance_tolerance = 0) {
  out <- .fit_climniche_matrix(
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
    preprocess = preprocess,
    preprocess_correlation = preprocess_correlation,
    preprocess_min_sd = preprocess_min_sd,
    global_mean = global_mean,
    global_sd = global_sd,
    tolerance = tolerance,
    tolerance_quantile = tolerance_quantile,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance
  )
  out$call <- match.call()
  out
}

#' Fit climniche to raster data
#'
#' @param current Raster* object of current climate layers.
#' @param future Raster* object of future climate layers with the same geometry
#'   and variables as `current`. Named layers are matched before fitting.
#' @param occupied Optional RasterLayer with binary or continuous occurrence,
#'   range, or SDM suitability values.
#' @param occupied_threshold Values at or below this threshold receive zero
#'   reference weight. Values above it keep their original value.
#' @param domain Optional RasterLayer limiting cells where exposure is analysed.
#' @param domain_threshold Values greater than this threshold define the domain.
#' @inheritParams fit_climniche
#'
#' @return An object of class `climniche_fit` with RasterLayer outputs stored in
#'   `x$rasters`.
#'
#' @details
#' `fit_climniche_raster()` is the Raster* workflow for users working with the
#' `raster` package. It extracts aligned cell values, applies `domain` and
#' `occupied` rasters, calls the matrix workflow, and writes the fitted
#' quantities back to RasterLayer outputs.
#' @export
fit_climniche_raster <- function(current, future, occupied = NULL,
                                 occupied_threshold = NULL, domain = NULL,
                                 domain_threshold = 0,
                                 cnfa = NULL, center = NULL,
                                 sensitivity = NULL, A = NULL,
                                 metric = c("diag", "factor"),
                                 boundary = 0.95, scale = TRUE,
                                 preprocess = TRUE,
                                 preprocess_correlation = 0.95,
                                 preprocess_min_sd = 1e-08,
                                 global_mean = NULL, global_sd = NULL,
                                 tolerance = NULL,
                                 tolerance_quantile = 0.10,
                                 boundary_exceedance_tolerance = 0) {
  out <- .fit_climniche_raster(current = current, future = future,
                        occupied = occupied,
                        occupied_threshold = occupied_threshold,
                        domain = domain,
                        domain_threshold = domain_threshold,
                        cnfa = cnfa,
                        center = center,
                        sensitivity = sensitivity,
                        A = A,
                        metric = metric,
                        boundary = boundary,
                        scale = scale,
                        preprocess = preprocess,
                        preprocess_correlation = preprocess_correlation,
                        preprocess_min_sd = preprocess_min_sd,
                        global_mean = global_mean,
                        global_sd = global_sd,
                        tolerance = tolerance,
                        tolerance_quantile = tolerance_quantile,
                        boundary_exceedance_tolerance =
                          boundary_exceedance_tolerance)
  out$call <- match.call()
  out
}

#' Fit climniche to terra raster data
#'
#' @param current terra SpatRaster of current climate layers.
#' @param future terra SpatRaster of future climate layers with the same
#'   geometry and variables as `current`. Named layers are matched before
#'   fitting.
#' @param occupied Optional one layer SpatRaster with binary or continuous
#'   occurrence, range, or SDM suitability values.
#' @param occupied_threshold Values at or below this threshold receive zero
#'   reference weight. Values above it keep their original value.
#' @param domain Optional one layer SpatRaster limiting cells where exposure is
#'   analysed.
#' @param domain_threshold Values greater than this threshold define the domain.
#' @inheritParams fit_climniche
#'
#' @return An object of class `climniche_fit` with SpatRaster outputs stored in
#'   `x$rasters`.
#'
#' @details
#' `fit_climniche_terra()` is the SpatRaster workflow for users working with
#' `terra`. It extracts aligned cell values, applies `domain` and `occupied`
#' rasters, calls the matrix workflow, and writes the fitted quantities back to
#' SpatRaster outputs.
#' @export
fit_climniche_terra <- function(current, future, occupied = NULL,
                                occupied_threshold = NULL, domain = NULL,
                                domain_threshold = 0,
                                cnfa = NULL, center = NULL,
                                sensitivity = NULL, A = NULL,
                                metric = c("diag", "factor"),
                                boundary = 0.95, scale = TRUE,
                                preprocess = TRUE,
                                preprocess_correlation = 0.95,
                                preprocess_min_sd = 1e-08,
                                global_mean = NULL, global_sd = NULL,
                                tolerance = NULL,
                                tolerance_quantile = 0.10,
                                boundary_exceedance_tolerance = 0) {
  out <- .fit_climniche_terra(current = current, future = future,
                       occupied = occupied,
                       occupied_threshold = occupied_threshold,
                       domain = domain,
                       domain_threshold = domain_threshold,
                       cnfa = cnfa,
                       center = center,
                       sensitivity = sensitivity,
                       A = A,
                       metric = metric,
                       boundary = boundary,
                       scale = scale,
                       preprocess = preprocess,
                       preprocess_correlation = preprocess_correlation,
                       preprocess_min_sd = preprocess_min_sd,
                       global_mean = global_mean,
                       global_sd = global_sd,
                       tolerance = tolerance,
                       tolerance_quantile = tolerance_quantile,
                       boundary_exceedance_tolerance =
                         boundary_exceedance_tolerance)
  out$call <- match.call()
  out
}

#' Plot a climniche map
#'
#' @param x A fitted climniche object with raster outputs, a RasterLayer, or a
#'   terra SpatRaster.
#' @param metric Quantity to plot. Accepted values are
#'   `"climate_change_amount"` for Climatic Displacement,
#'   `"niche_distance_change"` for Niche Distance Shift,
#'   `"climate_reconfiguration"` for Climatic Reconfiguration,
#'   `"niche_boundary_exceedance"` for Niche Boundary Exceedance, and
#'   `"change_alignment"` for the signed ratio used internally. Legacy aliases
#'   `"composition_change"` and `"outside_niche_exceedance"` still work.
#' @param occupied Optional current reference RasterLayer or terra SpatRaster
#'   used for masking or overlaying the current reference distribution.
#' @param occupied_only If TRUE, mask the plotted raster to current occurrence
#'   or suitability cells with positive reference weight.
#' @param occupied_threshold Threshold used when `occupied` contains binary or
#'   continuous values. Values above the threshold keep their original value
#'   when used as an overlay or mask.
#' @param title Optional plot title. Use `FALSE` to suppress it.
#' @param midpoint Midpoint for the Niche Distance Shift colour scale.
#' @param limits Optional two-element colour scale limits.
#' @param breaks Optional colour scale breaks.
#' @param colours Optional colour vector replacing the metric palette.
#' @param legend_title Optional colour legend title. The default, `FALSE`,
#'   suppresses a title that would repeat the panel title.
#' @param legend_position Position passed to the ggplot theme. Common values
#'   are `"right"`, `"bottom"`, and `"none"`.
#' @param show_legend If FALSE, suppress the colour legend.
#' @param symmetric If TRUE, use limits symmetric around `midpoint`.
#' @param extent Optional `c(xmin, xmax, ymin, ymax)` plotting extent.
#' @param degree_labels `"auto"` uses hemisphere degree labels for
#'   longitude-latitude rasters, `"hemisphere"` always uses them, and `"none"`
#'   uses the default ggplot labels.
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
#' @param colour_by Quantity used for point colour. The current plotting method
#'   uses `"niche_boundary_exceedance"`.
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
#' @param metric Quantity to plot. Accepted values are
#'   `"climate_change_amount"`, `"niche_distance_change"`,
#'   `"climate_reconfiguration"`, and `"niche_boundary_exceedance"`. Legacy
#'   aliases `"composition_change"` and `"outside_niche_exceedance"` still work.
#' @param scope `"current"` for a reference-weighted distribution or `"all"`
#'   for an unweighted distribution across evaluated cells.
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
