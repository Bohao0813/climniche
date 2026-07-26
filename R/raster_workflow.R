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
  current_complete <- rowSums(!is.finite(x0)) == 0L
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
    current_complete <- current_complete & is.finite(dom) &
      dom > domain_threshold
  }
  complete <- current_complete & rowSums(!is.finite(x1)) == 0L
  if (!any(complete)) {
    stop("No analysable cells remain after applying environmental and domain masks.",
         call. = FALSE)
  }

  if (is.null(occupied)) {
    occupied_values <- NULL
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
    occupied_values <- raster::getValues(occupied)
  }

  dots <- list(...)
  reference_fit <- do.call(
    .fit_climniche_matrix,
    c(
      list(
        current = x0[current_complete, , drop = FALSE],
        future = x0[current_complete, , drop = FALSE],
        occupied = if (is.null(occupied_values)) {
          NULL
        } else {
          occupied_values[current_complete]
        },
        occupied_threshold = occupied_threshold
      ),
      dots
    )
  )
  reference <- .reference_from_fit(reference_fit)
  descriptor_names <- c(
    "tolerance", "tolerance_quantile",
    "boundary_exceedance_tolerance"
  )
  projection_args <- c(
    list(
      reference = reference,
      current = x0[complete, , drop = FALSE],
      future = x1[complete, , drop = FALSE],
      occupied = if (is.null(occupied_values)) {
        NULL
      } else {
        occupied_values[complete]
      },
      occupied_threshold = occupied_threshold
    ),
    dots[intersect(names(dots), descriptor_names)]
  )
  fit <- do.call(project_climniche, projection_args)

  fit <- .attach_series_spatial_outputs(
    fit = fit,
    template = current,
    complete = complete,
    type = "raster"
  )
  fit$reference_cell_index <- which(current_complete)
  fit$reference_raster_complete <- current_complete
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
  current_complete <- rowSums(!is.finite(x0)) == 0L
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
    current_complete <- current_complete & is.finite(dom) &
      dom > domain_threshold
  }
  complete <- current_complete & rowSums(!is.finite(x1)) == 0L
  if (!any(complete)) {
    stop("No analysable cells remain after applying environmental and domain masks.",
         call. = FALSE)
  }

  if (is.null(occupied)) {
    occupied_values <- NULL
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
    occupied_values <- terra::values(occupied)[, 1]
  }

  dots <- list(...)
  reference_fit <- do.call(
    .fit_climniche_matrix,
    c(
      list(
        current = x0[current_complete, , drop = FALSE],
        future = x0[current_complete, , drop = FALSE],
        occupied = if (is.null(occupied_values)) {
          NULL
        } else {
          occupied_values[current_complete]
        },
        occupied_threshold = occupied_threshold
      ),
      dots
    )
  )
  reference <- .reference_from_fit(reference_fit)
  descriptor_names <- c(
    "tolerance", "tolerance_quantile",
    "boundary_exceedance_tolerance"
  )
  projection_args <- c(
    list(
      reference = reference,
      current = x0[complete, , drop = FALSE],
      future = x1[complete, , drop = FALSE],
      occupied = if (is.null(occupied_values)) {
        NULL
      } else {
        occupied_values[complete]
      },
      occupied_threshold = occupied_threshold
    ),
    dots[intersect(names(dots), descriptor_names)]
  )
  fit <- do.call(project_climniche, projection_args)

  fit <- .attach_series_spatial_outputs(
    fit = fit,
    template = current,
    complete = complete,
    type = "terra"
  )
  fit$reference_cell_index <- which(current_complete)
  fit$reference_raster_complete <- current_complete
  fit
}

.values_to_spatraster <- function(template, values_complete, complete) {
  r <- template[[1]]
  out <- rep(NA_real_, terra::ncell(r))
  out[complete] <- values_complete
  terra::values(r) <- out
  r
}
