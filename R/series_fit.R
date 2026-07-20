#' Fit niche relative climate exposure through time
#'
#' Fits one current climatic niche reference and evaluates an ordered set of
#' future projections. Projections may represent time periods, climate models,
#' scenarios, or their combinations.
#'
#' @param current Current climatic conditions supplied as a numeric matrix,
#'   data frame, Raster* object or SpatRaster.
#' @param future A list of future objects matching `current`. A three-dimensional
#'   numeric array may also be supplied, with cells, variables and projections
#'   in its three dimensions.
#' @param time Numeric, Date or POSIXct projection times, one per future object.
#' @param model Optional climate model identifiers, one per future object or a
#'   single value recycled across projections.
#' @param scenario Optional scenario identifiers, one per future object or a
#'   single value recycled across projections.
#' @param occupied Current reference rows, weights or a matching one-layer
#'   raster.
#' @param occupied_threshold Optional cutoff for continuous reference weights.
#' @param domain Optional one-layer raster limiting a spatial analysis.
#' @param domain_threshold Threshold used when `domain` is supplied.
#' @param cnfa Optional compatible CENFA object.
#' @param center Optional realised niche centre in the fitted climatic space.
#' @param sensitivity Optional non-negative variable sensitivity weights.
#' @param A Optional climatic weighting matrix.
#' @param metric Method used to build `A` when it is not supplied.
#' @param boundary Weighted quantile defining the empirical niche boundary.
#' @param scale If `TRUE`, standardise retained variables using current climate.
#' @param preprocess If `TRUE`, screen near-zero variance and highly correlated
#'   current-climate variables.
#' @param preprocess_correlation Maximum absolute correlation retained during
#'   preprocessing.
#' @param preprocess_min_sd Minimum current-climate standard deviation retained
#'   during preprocessing.
#' @param global_mean,global_sd Optional centring and scaling values.
#' @param tolerance Optional tolerance around zero for Niche Distance Shift.
#' @param tolerance_quantile Quantile used when `tolerance = NULL`.
#' @param boundary_exceedance_tolerance Tolerance used by the empirical boundary
#'   descriptor and dynamic summaries.
#'
#' @return A `climniche_series` object containing one compatible
#'   `climniche_fit` per projection.
#'
#' @details
#' The current niche centre, climatic weighting matrix, standardisation and
#' empirical boundary are held fixed. For spatial inputs, this reference is
#' estimated from finite current cells within `domain`, independently of
#' missing values in future projections. Future comparisons use the common set
#' of finite cells across all projections. All projections also use one Niche
#' Distance Shift tolerance. When `tolerance = NULL`, it is the requested
#' quantile of pooled absolute Niche Distance Shift values across the fitted
#' series.
#'
#' @examples
#' sim <- simulate_climniche(n = 180, p = 6, seed = 12)
#' future <- lapply(c(0.25, 0.50, 0.75, 1), function(fraction) {
#'   sim$current + fraction * (sim$future_away - sim$current)
#' })
#' series <- fit_climniche_series(
#'   current = sim$current,
#'   future = future,
#'   time = c(2030, 2050, 2070, 2090),
#'   occupied = sim$occupied,
#'   sensitivity = sim$sensitivity
#' )
#' climniche_range_summary(series)
#' @export
fit_climniche_series <- function(
    current,
    future,
    time,
    model = NULL,
    scenario = NULL,
    occupied = NULL,
    occupied_threshold = NULL,
    domain = NULL,
    domain_threshold = 0,
    cnfa = NULL,
    center = NULL,
    sensitivity = NULL,
    A = NULL,
    metric = c("diag", "factor"),
    boundary = 0.95,
    scale = TRUE,
    preprocess = TRUE,
    preprocess_correlation = 0.95,
    preprocess_min_sd = 1e-08,
    global_mean = NULL,
    global_sd = NULL,
    tolerance = NULL,
    tolerance_quantile = 0.10,
    boundary_exceedance_tolerance = 0) {
  metric <- match.arg(metric)
  projections <- .as_future_list(future, time)
  index <- .series_index(
    time = time,
    model = model,
    scenario = scenario,
    projection_names = names(projections)
  )
  names(projections) <- index$projection

  input_type <- .series_input_type(current)
  if (identical(input_type, "matrix")) {
    if (!is.null(domain)) {
      stop("domain is only used by raster and terra series.", call. = FALSE)
    }
    if (!all(vapply(projections, .is_matrix_input, logical(1)))) {
      stop("All future projections must be numeric matrices or data frames.",
           call. = FALSE)
    }
    reference <- fit_climniche_reference(
      current = current,
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
      global_sd = global_sd
    )
    fits <- lapply(projections, function(future_projection) {
      project_climniche(
        reference = reference,
        future = future_projection,
        current = current,
        occupied = occupied,
        occupied_threshold = occupied_threshold,
        tolerance = tolerance,
        tolerance_quantile = tolerance_quantile,
        boundary_exceedance_tolerance =
          boundary_exceedance_tolerance
      )
    })
  } else {
    .check_spatial_series_types(projections, input_type)
    fit_fun <- if (identical(input_type, "raster")) {
      fit_climniche_raster
    } else {
      fit_climniche_terra
    }
    reference_fit <- do.call(fit_fun, list(
      current = current,
      future = current,
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
      tolerance = 0,
      tolerance_quantile = tolerance_quantile,
      boundary_exceedance_tolerance = 0
    ))
    reference <- .reference_from_fit(reference_fit)
    common_domain <- .series_common_domain(
      current = current,
      projections = projections,
      domain = domain,
      domain_threshold = domain_threshold,
      type = input_type
    )
    complete <- .series_domain_cells(common_domain, input_type)
    current_values <- .series_spatial_values(current, input_type)
    reference_weight <- rep(0, nrow(current_values))
    reference_weight[reference_fit$cell_index] <- reference$reference_weight
    projected_weight <- reference_weight[complete]
    fits <- lapply(projections, function(future_projection) {
      future_values <- .series_spatial_values(future_projection, input_type)
      fit <- project_climniche(
        reference = reference,
        current = current_values[complete, , drop = FALSE],
        future = future_values[complete, , drop = FALSE],
        occupied = projected_weight,
        tolerance = tolerance,
        tolerance_quantile = tolerance_quantile,
        boundary_exceedance_tolerance = boundary_exceedance_tolerance
      )
      .attach_series_spatial_outputs(
        fit = fit,
        template = current,
        complete = complete,
        type = input_type
      )
    })
  }
  fits <- .harmonize_series_descriptors(
    fits = fits,
    tolerance = tolerance,
    tolerance_quantile = tolerance_quantile,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance
  )
  names(fits) <- index$projection
  .validate_series_fits(fits)

  for (i in seq_along(fits)) {
    fits[[i]]$projection <- index$projection[i]
    fits[[i]]$projection_time <- index$time[i]
    fits[[i]]$climate_model <- index$model[i]
    fits[[i]]$climate_scenario <- index$scenario[i]
  }
  out <- list(
    call = match.call(),
    fits = fits,
    index = index,
    reference = reference,
    descriptor_settings = fits[[1L]]$descriptor_settings,
    input_type = input_type
  )
  class(out) <- "climniche_series"
  out
}

.series_spatial_values <- function(x, type) {
  if (identical(type, "raster")) {
    return(.raster_values_matrix(x))
  }
  terra::values(x, mat = TRUE)
}

.series_domain_cells <- function(x, type) {
  values <- if (identical(type, "raster")) {
    raster::getValues(x)
  } else {
    terra::values(x)[, 1L]
  }
  is.finite(values) & values > 0
}

.attach_series_spatial_outputs <- function(fit, template, complete, type) {
  converter <- if (identical(type, "raster")) {
    .values_to_raster
  } else {
    .values_to_spatraster
  }
  metric_values <- list(
    climate_change_amount = fit$climate_change_amount,
    niche_distance_change = fit$niche_distance_change,
    climate_reconfiguration = fit$climate_reconfiguration,
    composition_change = fit$composition_change,
    change_alignment = fit$change_alignment,
    niche_boundary_exceedance = fit$niche_boundary_exceedance,
    outside_niche_exceedance = fit$outside_niche_exceedance,
    radial_direction = as.integer(fit$radial_direction),
    boundary_status = as.integer(fit$boundary_status)
  )
  fit$rasters <- lapply(metric_values, function(values) {
    converter(template, values, complete)
  })
  for (name in names(fit$rasters)) {
    names(fit$rasters[[name]]) <- name
  }
  fit$descriptor_lookup <- list(
    radial_direction = data.frame(
      id = seq_along(levels(fit$radial_direction)),
      level = levels(fit$radial_direction)
    ),
    boundary_status = data.frame(
      id = seq_along(levels(fit$boundary_status)),
      level = levels(fit$boundary_status)
    )
  )
  fit$cell_index <- which(complete)
  fit$raster_complete <- complete
  fit
}

.harmonize_series_descriptors <- function(
    fits,
    tolerance,
    tolerance_quantile,
    boundary_exceedance_tolerance) {
  tolerance_quantile <- .check_probability(
    tolerance_quantile,
    "tolerance_quantile"
  )
  if (is.null(tolerance)) {
    distance_change <- unlist(lapply(
      fits,
      function(fit) abs(fit$niche_distance_change)
    ), use.names = FALSE)
    tolerance <- as.numeric(stats::quantile(
      distance_change,
      probs = tolerance_quantile,
      names = FALSE,
      na.rm = TRUE,
      type = 8
    ))
  } else {
    tolerance <- .check_finite_scalar(tolerance, "tolerance")
  }
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a finite non-negative number.", call. = FALSE)
  }

  lapply(fits, function(fit) {
    descriptors <- .exposure_descriptors(
      niche_distance_change = fit$niche_distance_change,
      niche_boundary_exceedance = fit$niche_boundary_exceedance,
      tolerance = tolerance,
      tolerance_quantile = tolerance_quantile,
      boundary_exceedance_tolerance = boundary_exceedance_tolerance
    )
    settings <- attr(descriptors, "descriptor_settings")
    attr(descriptors, "descriptor_settings") <- NULL
    fit$radial_direction <- descriptors$radial_direction
    fit$boundary_status <- descriptors$boundary_status
    fit$descriptor_settings <- settings
    fit$threshold_settings <- settings

    if (!is.null(fit$rasters) && !is.null(fit$raster_complete)) {
      template <- .fit_spatial_template(fit)
      converter <- if (.is_raster(template)) {
        .values_to_raster
      } else {
        .values_to_spatraster
      }
      fit$rasters$radial_direction <- converter(
        template,
        as.integer(fit$radial_direction),
        fit$raster_complete
      )
      fit$rasters$boundary_status <- converter(
        template,
        as.integer(fit$boundary_status),
        fit$raster_complete
      )
      names(fit$rasters$radial_direction) <- "radial_direction"
      names(fit$rasters$boundary_status) <- "boundary_status"
    }
    fit
  })
}

.is_matrix_input <- function(x) {
  is.matrix(x) || is.data.frame(x)
}

.series_input_type <- function(current) {
  if (.is_matrix_input(current)) {
    return("matrix")
  }
  if (.is_raster(current)) {
    return("raster")
  }
  if (.is_spatraster(current)) {
    return("terra")
  }
  stop("current must be a matrix, data frame, Raster* object or SpatRaster.",
       call. = FALSE)
}

.as_future_list <- function(future, time) {
  n_time <- length(time)
  if (is.array(future) && length(dim(future)) == 3L) {
    projections <- lapply(seq_len(dim(future)[3L]), function(i) {
      out <- future[, , i, drop = FALSE]
      dim(out) <- dim(future)[1:2]
      dimnames(out) <- dimnames(future)[1:2]
      out
    })
    names(projections) <- dimnames(future)[[3L]]
  } else if (is.list(future) && !is.data.frame(future)) {
    projections <- future
  } else if (n_time == 1L) {
    projections <- list(future)
  } else {
    stop("future must be a list or a three-dimensional numeric array.",
         call. = FALSE)
  }
  if (!length(projections) || length(projections) != n_time) {
    stop("time must contain one value per future projection.", call. = FALSE)
  }
  projections
}

.series_index <- function(time, model, scenario, projection_names = NULL) {
  if (!(is.numeric(time) || inherits(time, "Date") ||
        inherits(time, "POSIXct") || inherits(time, "POSIXlt"))) {
    stop("time must be numeric, Date or POSIXct.", call. = FALSE)
  }
  if (any(is.na(time)) || any(!is.finite(as.numeric(time)))) {
    stop("time must contain finite, non-missing values.", call. = FALSE)
  }
  n <- length(time)
  model <- .recycle_series_label(model, n, "model")
  scenario <- .recycle_series_label(scenario, n, "scenario")
  key <- paste(
    as.numeric(time),
    .series_label_key(model),
    .series_label_key(scenario),
    sep = "\r"
  )
  if (anyDuplicated(key)) {
    stop("Each time, model and scenario combination must be unique.",
         call. = FALSE)
  }
  if (is.null(projection_names) ||
      length(projection_names) != n ||
      anyNA(projection_names) ||
      any(!nzchar(projection_names)) ||
      anyDuplicated(projection_names)) {
    projection_names <- sprintf("projection_%03d", seq_len(n))
  }
  data.frame(
    projection = projection_names,
    time = time,
    model = model,
    scenario = scenario,
    stringsAsFactors = FALSE
  )
}

.recycle_series_label <- function(x, n, arg) {
  if (is.null(x)) {
    return(rep(NA_character_, n))
  }
  if (length(x) == 1L) {
    x <- rep(x, n)
  }
  if (length(x) != n) {
    stop(arg, " must have length one or one value per projection.",
         call. = FALSE)
  }
  x <- as.character(x)
  x[!nzchar(x)] <- NA_character_
  x
}

.series_label_key <- function(x) {
  ifelse(is.na(x), "<unspecified>", x)
}

.check_spatial_series_types <- function(projections, type) {
  ok <- if (identical(type, "raster")) {
    vapply(projections, .is_raster, logical(1))
  } else {
    vapply(projections, .is_spatraster, logical(1))
  }
  if (!all(ok)) {
    stop("All future projections must use the same spatial class as current.",
         call. = FALSE)
  }
}

.series_common_domain <- function(current, projections, domain,
                                  domain_threshold, type) {
  domain_threshold <- .check_finite_scalar(
    domain_threshold,
    "domain_threshold"
  )
  if (identical(type, "raster")) {
    .need_raster()
    current_values <- .raster_values_matrix(current)
    complete <- rowSums(!is.finite(current_values)) == 0L
    for (future in projections) {
      if (raster::nlayers(current) != raster::nlayers(future) ||
          !raster::compareRaster(current, future, stopiffalse = FALSE)) {
        stop("All future rasters must match current geometry and layers.",
             call. = FALSE)
      }
      values <- .raster_values_matrix(future)
      complete <- complete & rowSums(!is.finite(values)) == 0L
    }
    if (!is.null(domain)) {
      if (!.is_raster(domain) || raster::nlayers(domain) != 1L ||
          !raster::compareRaster(raster::raster(current), domain,
                                 stopiffalse = FALSE)) {
        stop("domain must be a one-layer raster matching current.",
             call. = FALSE)
      }
      values <- raster::getValues(domain)
      complete <- complete & is.finite(values) & values > domain_threshold
    }
    if (!any(complete)) {
      stop("No common analysable raster cells remain across projections.",
           call. = FALSE)
    }
    out <- raster::raster(current)
    raster::values(out) <- as.numeric(complete)
    return(out)
  }

  .need_terra()
  current_values <- terra::values(current, mat = TRUE)
  complete <- rowSums(!is.finite(current_values)) == 0L
  for (future in projections) {
    if (terra::nlyr(current) != terra::nlyr(future) ||
        !terra::compareGeom(current, future, stopOnError = FALSE)) {
      stop("All future rasters must match current geometry and layers.",
           call. = FALSE)
    }
    values <- terra::values(future, mat = TRUE)
    complete <- complete & rowSums(!is.finite(values)) == 0L
  }
  if (!is.null(domain)) {
    if (!.is_spatraster(domain) || terra::nlyr(domain) != 1L ||
        !terra::compareGeom(current[[1L]], domain, stopOnError = FALSE)) {
      stop("domain must be a one-layer SpatRaster matching current.",
           call. = FALSE)
    }
    values <- terra::values(domain)[, 1L]
    complete <- complete & is.finite(values) & values > domain_threshold
  }
  if (!any(complete)) {
    stop("No common analysable raster cells remain across projections.",
         call. = FALSE)
  }
  out <- current[[1L]]
  terra::values(out) <- as.numeric(complete)
  out
}

.validate_series_fits <- function(fits) {
  if (!length(fits) || !all(vapply(
    fits,
    inherits,
    logical(1),
    what = "climniche_fit"
  ))) {
    stop("Series projections did not produce compatible climniche fits.",
         call. = FALSE)
  }
  first <- fits[[1L]]
  for (i in seq_along(fits)[-1L]) {
    candidate <- fits[[i]]
    same <- isTRUE(all.equal(first$center, candidate$center, tolerance = 0)) &&
      isTRUE(all.equal(first$A, candidate$A, tolerance = 0)) &&
      isTRUE(all.equal(first$boundary_radius,
                       candidate$boundary_radius, tolerance = 0)) &&
      identical(colnames(first$current), colnames(candidate$current)) &&
      identical(first$cell_index, candidate$cell_index)
    if (!same) {
      stop("All series projections must use one fixed climatic niche reference and cell support.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Extract a long table from a climniche series
#'
#' @param x A `climniche_series` object.
#' @param scope `"current"` for positive reference weights or `"all"` for all
#'   evaluated cells.
#'
#' @return A data frame with one row per cell and projection.
#' @export
climniche_series_table <- function(x, scope = c("current", "all")) {
  if (!inherits(x, "climniche_series")) {
    stop("x must be a climniche_series object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  rows <- lapply(seq_along(x$fits), function(i) {
    tab <- climniche_table(x$fits[[i]], scope = scope)
    meta <- x$index[i, , drop = FALSE]
    data.frame(
      projection = rep(meta$projection, nrow(tab)),
      time = rep(meta$time, nrow(tab)),
      model = rep(meta$model, nrow(tab)),
      scenario = rep(meta$scenario, nrow(tab)),
      tab,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.series_range_summary <- function(
    x,
    scope = c("current", "all"),
    aggregation_weight = NULL,
    area_weight = FALSE,
    boundary_exceedance_tolerance = NULL) {
  scope <- match.arg(scope)
  rows <- lapply(seq_along(x$fits), function(i) {
    summary <- .fit_range_summary(
      x$fits[[i]],
      scope = scope,
      aggregation_weight = aggregation_weight,
      area_weight = area_weight,
      boundary_exceedance_tolerance =
        boundary_exceedance_tolerance
    )
    meta <- x$index[i, , drop = FALSE]
    data.frame(meta, summary, check.names = FALSE,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @export
print.climniche_series <- function(x, ...) {
  cat("Niche relative climate exposure series\n")
  cat("Projections:", nrow(x$index), "\n")
  cat("Times:", length(unique(x$index$time)), "\n")
  model_count <- length(unique(x$index$model[!is.na(x$index$model)]))
  scenario_count <- length(unique(
    x$index$scenario[!is.na(x$index$scenario)]
  ))
  cat("Climate models:", if (model_count) model_count else 1L, "\n")
  cat("Scenarios:", if (scenario_count) scenario_count else 1L, "\n")
  cat("Evaluated cells:", nrow(x$fits[[1L]]$current), "\n")
  invisible(x)
}
