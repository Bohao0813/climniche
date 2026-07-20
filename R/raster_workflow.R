#' Fit climniche to Raster* climate layers
#'
#' @param current Raster* object of current climate layers.
#' @param future Raster* object of future climate layers with the same geometry
#'   and variables as `current`. Named layers are matched before fitting.
#' @param occupied Optional RasterLayer indicating current occurrence or range
#'   cells, or continuous SDM suitability values used as reference weights.
#' @param occupied_threshold Numeric threshold used when `occupied` is a raster.
#' @param domain Optional RasterLayer limiting the cells to analyse. Values
#'   greater than `domain_threshold` define the analysis domain.
#' @param domain_threshold Numeric threshold used when `domain` is supplied.
#' @param ... Arguments passed to `fit_climniche()`.
#'
#' @return A `climniche_fit` object with an additional `rasters` element containing
#'   RasterLayer outputs.
#' @noRd
.fit_climniche_raster <- function(current, future, occupied = NULL,
                                  occupied_threshold = NULL, domain = NULL,
                                  domain_threshold = 0, ...) {
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("The raster package is required for fit_climniche_raster().",
         call. = FALSE)
  }
  if (!methods::is(current, "Raster") || !methods::is(future, "Raster")) {
    stop("current and future must be raster Raster* objects.", call. = FALSE)
  }
  domain_threshold <- .check_finite_scalar(domain_threshold,
                                            "domain_threshold")
  if (raster::nlayers(current) < 1L ||
      raster::nlayers(current) != raster::nlayers(future)) {
    stop("current and future rasters must have the same positive number of layers.",
         call. = FALSE)
  }
  if (!raster::compareRaster(current, future, stopiffalse = FALSE)) {
    stop("current and future rasters must have matching geometry.", call. = FALSE)
  }

  x0 <- .raster_values_matrix(current)
  x1 <- .raster_values_matrix(future)
  complete <- rowSums(!is.finite(x0)) == 0L &
    rowSums(!is.finite(x1)) == 0L
  if (!is.null(domain)) {
    if (!methods::is(domain, "Raster")) {
      stop("domain must be NULL or a RasterLayer.", call. = FALSE)
    }
    if (raster::nlayers(domain) != 1L) {
      stop("domain must have one layer.", call. = FALSE)
    }
    if (!raster::compareRaster(raster::raster(current), domain,
                               stopiffalse = FALSE)) {
      stop("domain raster must match current raster geometry.",
           call. = FALSE)
    }
    dom <- raster::getValues(domain)
    complete <- complete & is.finite(dom) & dom > domain_threshold
  }
  if (!any(complete)) {
    stop("No analysable cells remain after applying environmental and domain masks.",
         call. = FALSE)
  }

  if (is.null(occupied)) {
    occupied_weight <- NULL
  } else {
    if (!methods::is(occupied, "Raster")) {
      stop("occupied must be NULL or a RasterLayer.", call. = FALSE)
    }
    if (raster::nlayers(occupied) != 1L) {
      stop("occupied must have one layer.", call. = FALSE)
    }
    if (!raster::compareRaster(raster::raster(current), occupied,
                               stopiffalse = FALSE)) {
      stop("occupied raster must match current raster geometry.", call. = FALSE)
    }
    occ <- raster::getValues(occupied)
    occupied_weight <- occ[complete]
  }

  fit <- .fit_climniche_matrix(
    current = x0[complete, , drop = FALSE],
    future = x1[complete, , drop = FALSE],
    occupied = occupied_weight,
    occupied_threshold = occupied_threshold,
    ...
  )

  fit$rasters <- list(
    climate_change_amount = .values_to_raster(current, fit$climate_change_amount,
                                              complete),
    niche_distance_change = .values_to_raster(current, fit$niche_distance_change,
                                             complete),
    climate_reconfiguration = .values_to_raster(
      current,
      fit$climate_reconfiguration,
      complete
    ),
    composition_change = .values_to_raster(
      current,
      fit$composition_change,
      complete
    ),
    change_alignment = .values_to_raster(
      current,
      fit$change_alignment,
      complete
    ),
    niche_boundary_exceedance = .values_to_raster(
      current,
      fit$niche_boundary_exceedance,
      complete
    ),
    outside_niche_exceedance = .values_to_raster(
      current,
      fit$outside_niche_exceedance,
      complete
    ),
    radial_direction = .values_to_raster(
      current,
      as.integer(fit$radial_direction),
      complete
    ),
    boundary_status = .values_to_raster(
      current,
      as.integer(fit$boundary_status),
      complete
    )
  )
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
  for (nm in names(fit$rasters)) {
    names(fit$rasters[[nm]]) <- nm
  }
  fit$cell_index <- which(complete)
  fit$raster_complete <- complete
  fit
}

.raster_values_matrix <- function(x) {
  values <- raster::getValues(x)
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1L)
    colnames(values) <- names(x)
  }
  values
}

.values_to_raster <- function(template, values_complete, complete) {
  r <- raster::raster(template)
  out <- rep(NA_real_, raster::ncell(r))
  out[complete] <- values_complete
  raster::values(r) <- out
  r
}

#' Fit climniche to terra SpatRaster environmental layers
#'
#' @param current SpatRaster of current environmental layers.
#' @param future SpatRaster of future environmental layers with the same
#'   geometry and variables as `current`. Named layers are matched before
#'   fitting.
#' @param occupied Optional SpatRaster layer indicating current occurrence,
#'   range, suitability or probability values used as reference weights.
#' @param occupied_threshold Numeric threshold used when `occupied` is supplied.
#' @param domain Optional one-layer SpatRaster limiting the cells to analyse.
#'   Values greater than `domain_threshold` define the analysis domain.
#' @param domain_threshold Numeric threshold used when `domain` is supplied.
#' @param ... Arguments passed to `fit_climniche()`.
#'
#' @return A `climniche_fit` object with an additional `rasters` element containing
#'   SpatRaster outputs.
#' @noRd
.fit_climniche_terra <- function(current, future, occupied = NULL,
                                 occupied_threshold = NULL, domain = NULL,
                                 domain_threshold = 0, ...) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("The terra package is required for fit_climniche_terra().",
         call. = FALSE)
  }
  if (!methods::is(current, "SpatRaster") ||
      !methods::is(future, "SpatRaster")) {
    stop("current and future must be terra SpatRaster objects.",
         call. = FALSE)
  }
  domain_threshold <- .check_finite_scalar(domain_threshold,
                                            "domain_threshold")
  if (terra::nlyr(current) < 1L ||
      terra::nlyr(current) != terra::nlyr(future)) {
    stop("current and future rasters must have the same positive number of layers.",
         call. = FALSE)
  }
  if (!terra::compareGeom(current, future, stopOnError = FALSE)) {
    stop("current and future rasters must have matching geometry.",
         call. = FALSE)
  }

  x0 <- terra::values(current, mat = TRUE)
  x1 <- terra::values(future, mat = TRUE)
  complete <- rowSums(!is.finite(x0)) == 0L &
    rowSums(!is.finite(x1)) == 0L
  if (!is.null(domain)) {
    if (!methods::is(domain, "SpatRaster")) {
      stop("domain must be NULL or a terra SpatRaster.", call. = FALSE)
    }
    if (terra::nlyr(domain) != 1) {
      stop("domain must have one layer.", call. = FALSE)
    }
    if (!terra::compareGeom(current[[1]], domain, stopOnError = FALSE)) {
      stop("domain raster must match current raster geometry.",
           call. = FALSE)
    }
    dom <- terra::values(domain)[, 1]
    complete <- complete & is.finite(dom) & dom > domain_threshold
  }
  if (!any(complete)) {
    stop("No analysable cells remain after applying environmental and domain masks.",
         call. = FALSE)
  }

  if (is.null(occupied)) {
    occupied_weight <- NULL
  } else {
    if (!methods::is(occupied, "SpatRaster")) {
      stop("occupied must be NULL or a terra SpatRaster.", call. = FALSE)
    }
    if (terra::nlyr(occupied) != 1) {
      stop("occupied must have one layer.", call. = FALSE)
    }
    if (!terra::compareGeom(current[[1]], occupied, stopOnError = FALSE)) {
      stop("occupied raster must match current raster geometry.",
           call. = FALSE)
    }
    occ <- terra::values(occupied)[, 1]
    occupied_weight <- occ[complete]
  }

  fit <- .fit_climniche_matrix(
    current = x0[complete, , drop = FALSE],
    future = x1[complete, , drop = FALSE],
    occupied = occupied_weight,
    occupied_threshold = occupied_threshold,
    ...
  )

  fit$rasters <- list(
    climate_change_amount = .values_to_spatraster(
      current, fit$climate_change_amount, complete
    ),
    niche_distance_change = .values_to_spatraster(
      current, fit$niche_distance_change, complete
    ),
    climate_reconfiguration = .values_to_spatraster(
      current, fit$climate_reconfiguration, complete
    ),
    composition_change = .values_to_spatraster(
      current, fit$composition_change, complete
    ),
    change_alignment = .values_to_spatraster(
      current, fit$change_alignment, complete
    ),
    niche_boundary_exceedance = .values_to_spatraster(
      current, fit$niche_boundary_exceedance, complete
    ),
    outside_niche_exceedance = .values_to_spatraster(
      current, fit$outside_niche_exceedance, complete
    ),
    radial_direction = .values_to_spatraster(
      current, as.integer(fit$radial_direction), complete
    ),
    boundary_status = .values_to_spatraster(
      current, as.integer(fit$boundary_status), complete
    )
  )
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
  for (nm in names(fit$rasters)) {
    names(fit$rasters[[nm]]) <- nm
  }
  fit$cell_index <- which(complete)
  fit$raster_complete <- complete
  fit
}

.values_to_spatraster <- function(template, values_complete, complete) {
  r <- template[[1]]
  out <- rep(NA_real_, terra::ncell(r))
  out[complete] <- values_complete
  terra::values(r) <- out
  r
}
