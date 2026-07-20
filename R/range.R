#' Summarise niche boundary exceedance across a range
#'
#' @param x A `climniche_fit` or `climniche_series` object.
#' @param scope `"current"` restricts the summary to positive reference
#'   weights. `"all"` uses every evaluated cell.
#' @param aggregation_weight Optional non-negative cell weights used only for
#'   the range summary. Numeric vectors and, for spatial fits, matching
#'   RasterLayer or SpatRaster objects are accepted. These weights multiply the
#'   reference weights when `scope = "current"`.
#' @param area_weight If `TRUE`, multiply the summary weights by raster cell
#'   area. This option requires a spatial fit.
#' @param boundary_exceedance_tolerance Optional non-negative tolerance used to
#'   identify positive Niche Boundary Exceedance. By default, the fitted
#'   descriptor tolerance is used.
#'
#' @return A data frame. For a `climniche_series`, one row is returned for each
#'   time, model and scenario combination.
#'
#' @details
#' Let `a_i` be the range summary weight, let \eqn{\tau} be the boundary
#' tolerance, and define
#' \deqn{\widetilde{E}_i = E_i I(E_i > \tau), \qquad
#'       e_i = \widetilde{E}_i / B_q.}
#' The weighted exposed fraction is
#' \deqn{F = \frac{\sum_i a_i I(e_i > 0)}{\sum_i a_i}.}
#' Conditional relative exceedance is
#' \deqn{S = \frac{\sum_i a_i e_i I(e_i > 0)}
#'                  {\sum_i a_i I(e_i > 0)}.}
#' The range mean relative exceedance is
#' \deqn{X = \frac{\sum_i a_i e_i I(e_i > 0)}{\sum_i a_i} = F S.}
#' For `scope = "current"`, `a_i` contains the reference weight. Optional
#' aggregation and cell area weights multiply it. For `scope = "all"`, the
#' reference weight is omitted. The result describes climatic exposure; it is
#' not an estimate of extinction risk or adaptive capacity.
#'
#' @examples
#' sim <- simulate_climniche(n = 250, p = 6, seed = 8)
#' fit <- fit_climniche(
#'   sim$current,
#'   sim$future_away,
#'   occupied = sim$occupied,
#'   sensitivity = sim$sensitivity
#' )
#' climniche_range_summary(fit)
#' @export
climniche_range_summary <- function(
    x,
    scope = c("current", "all"),
    aggregation_weight = NULL,
    area_weight = FALSE,
    boundary_exceedance_tolerance = NULL) {
  if (inherits(x, "climniche_series")) {
    return(.series_range_summary(
      x = x,
      scope = scope,
      aggregation_weight = aggregation_weight,
      area_weight = area_weight,
      boundary_exceedance_tolerance = boundary_exceedance_tolerance
    ))
  }
  .fit_range_summary(
    x = x,
    scope = scope,
    aggregation_weight = aggregation_weight,
    area_weight = area_weight,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance
  )
}

.fit_range_summary <- function(
    x,
    scope = c("current", "all"),
    aggregation_weight = NULL,
    area_weight = FALSE,
    boundary_exceedance_tolerance = NULL) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a climniche_fit or climniche_series object.",
         call. = FALSE)
  }
  scope <- match.arg(scope)
  area_weight <- .check_flag(area_weight, "area_weight")
  tolerance <- .range_boundary_tolerance(
    x,
    boundary_exceedance_tolerance
  )
  reference_weight <- .fit_reference_weights(x)
  modifier <- .aggregation_modifier(x, aggregation_weight)
  if (area_weight) {
    modifier <- modifier * .spatial_area_weights(x)
  }
  weights <- if (identical(scope, "current")) {
    reference_weight * modifier
  } else {
    modifier
  }
  weights <- .clean_aggregation_weights(weights, "summary weights")
  if (!any(weights > 0)) {
    stop("No positive range-summary weights remain.", call. = FALSE)
  }

  exceedance <- .fit_metric(x, "niche_boundary_exceedance")
  boundary <- as.numeric(x$boundary_radius)[1L]
  if (!is.finite(boundary) || boundary < 0) {
    stop("x contains an invalid niche boundary distance.", call. = FALSE)
  }
  exposed <- is.finite(exceedance) & exceedance > tolerance
  effective <- ifelse(exposed, exceedance, 0)
  total_weight <- sum(weights)
  exposed_weight <- sum(weights[exposed])
  exposed_fraction <- exposed_weight / total_weight
  mean_exceedance <- sum(weights * effective) / total_weight

  if (boundary > 0) {
    relative <- effective / boundary
    conditional <- if (exposed_weight > 0) {
      sum(weights * relative) / exposed_weight
    } else {
      0
    }
    range_wide <- sum(weights * relative) / total_weight
  } else {
    conditional <- NA_real_
    range_wide <- NA_real_
  }

  data.frame(
    scope = scope,
    n = sum(weights > 0),
    aggregation_weight_sum = total_weight,
    boundary_distance = boundary,
    boundary_exceedance_tolerance = tolerance,
    exposed_fraction = exposed_fraction,
    conditional_relative_exceedance = conditional,
    range_wide_relative_exceedance = range_wide,
    mean_niche_boundary_exceedance = mean_exceedance,
    stringsAsFactors = FALSE
  )
}

.range_boundary_tolerance <- function(x, tolerance) {
  if (is.null(tolerance)) {
    settings <- x$descriptor_settings %||% x$threshold_settings
    tolerance <- settings$boundary_exceedance_tolerance %||% 0
  }
  tolerance <- .check_finite_scalar(
    tolerance,
    "boundary_exceedance_tolerance"
  )
  if (tolerance < 0) {
    stop("boundary_exceedance_tolerance must be non-negative.",
         call. = FALSE)
  }
  tolerance
}

.clean_aggregation_weights <- function(x, arg = "aggregation_weight") {
  if (!is.numeric(x)) {
    stop(arg, " must be numeric.", call. = FALSE)
  }
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  if (any(x < 0)) {
    stop(arg, " contains negative values.", call. = FALSE)
  }
  x
}

.aggregation_modifier <- function(x, aggregation_weight) {
  n <- nrow(x$current)
  if (is.null(aggregation_weight)) {
    return(rep(1, n))
  }
  if (is.numeric(aggregation_weight)) {
    values <- aggregation_weight
    if (length(values) == n) {
      values <- .align_reference_values(
        values,
        rownames(x$current),
        arg = "aggregation_weight"
      )
      return(.clean_aggregation_weights(values))
    }
    if (!is.null(x$raster_complete) &&
        length(values) == length(x$raster_complete)) {
      return(.clean_aggregation_weights(values[x$raster_complete]))
    }
    stop("aggregation_weight must match the evaluated rows or raster cells.",
         call. = FALSE)
  }

  template <- .fit_spatial_template(x)
  if (is.null(template)) {
    stop("Raster aggregation weights require a spatial climniche fit.",
         call. = FALSE)
  }
  if (.is_raster(aggregation_weight)) {
    .need_raster()
    if (!.is_raster(template) || raster::nlayers(aggregation_weight) != 1L ||
        !raster::compareRaster(template, aggregation_weight,
                               stopiffalse = FALSE)) {
      stop("aggregation_weight raster must match the fitted raster geometry.",
           call. = FALSE)
    }
    values <- raster::getValues(aggregation_weight)
  } else if (.is_spatraster(aggregation_weight)) {
    .need_terra()
    if (!.is_spatraster(template) || terra::nlyr(aggregation_weight) != 1L ||
        !terra::compareGeom(template, aggregation_weight,
                            stopOnError = FALSE)) {
      stop("aggregation_weight raster must match the fitted raster geometry.",
           call. = FALSE)
    }
    values <- terra::values(aggregation_weight)[, 1L]
  } else {
    stop("aggregation_weight must be numeric or a matching one-layer raster.",
         call. = FALSE)
  }
  .clean_aggregation_weights(values[x$raster_complete])
}

.fit_spatial_template <- function(x) {
  if (is.null(x$rasters) || !length(x$rasters)) {
    return(NULL)
  }
  x$rasters[[1L]]
}

.spatial_area_weights <- function(x) {
  template <- .fit_spatial_template(x)
  if (is.null(template) || is.null(x$raster_complete)) {
    stop("area_weight = TRUE requires a spatial climniche fit.",
         call. = FALSE)
  }
  if (.is_raster(template)) {
    .need_raster()
    values <- raster::getValues(raster::area(template))
  } else if (.is_spatraster(template)) {
    .need_terra()
    values <- terra::values(terra::cellSize(template, unit = "km"))[, 1L]
  } else {
    stop("Unsupported spatial raster type.", call. = FALSE)
  }
  values <- values[x$raster_complete]
  values <- .clean_aggregation_weights(values, "cell areas")
  if (!any(values > 0)) {
    stop("No positive raster cell areas are available.", call. = FALSE)
  }
  values
}
