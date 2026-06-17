#' Fit niche relative climate exposure
#'
#' @param current Numeric matrix or data frame of current climate values.
#' @param future Numeric matrix or data frame of future climate values, with the
#'   same dimensions and variables as `current`.
#' @param occupied NULL, logical vector, row indices, or a numeric vector with
#'   one value per row. Numeric vectors of length `nrow(current)` are treated as
#'   continuous reference weights, including SDM suitability values.
#' @param occupied_threshold Threshold used when `occupied` is a numeric vector
#'   with one value per row.
#' @param cnfa Optional CENFA `cnfa` object. When supplied, its `mf` and `sf`
#'   slots are used unless `center`, `sensitivity`, or `A` are provided.
#' @param center Optional realised niche centre in standardized climate space.
#' @param sensitivity Optional climate-variable sensitivity weights.
#' @param A Optional niche metric matrix.
#' @param metric `"diag"` or `"factor"`. Used only when `A` is missing.
#' @param boundary Quantile for the empirical occupied niche boundary.
#' @param scale Logical. If TRUE, standardize current and future values using
#'   means and standard deviations of current values.
#' @param global_mean Optional means used for standardization.
#' @param global_sd Optional standard deviations used for standardization.
#' @param tolerance Optional Niche Distance Shift tolerance.
#' @param tolerance_quantile Quantile of absolute Niche Distance Shift used
#'   when `tolerance = NULL`.
#' @param stable_climate_change Optional Climatic Displacement threshold for
#'   the limited niche-relative change class.
#' @param stable_quantile Quantile of Climatic Displacement used when
#'   `stable_climate_change = NULL`.
#' @param stable_reconfiguration Optional threshold for low Climatic
#'   Reconfiguration.
#' @param stable_reconfiguration_quantile Quantile of Climatic Reconfiguration
#'   used when `stable_reconfiguration = NULL`.
#' @param boundary_exceedance_tolerance Tolerance for deciding whether future
#'   climate exceeds the empirical niche boundary.
#' @param conflict_ratio Minimum minority-sign variable-contribution share used
#'   to mark mixed variable responses. Set to NULL to disable this flag.
#'
#' @return An object of class `climniche_fit`.
#' @noRd
.fit_climniche_matrix <- function(current, future, occupied = NULL,
                                  occupied_threshold = NULL,
                                  cnfa = NULL, center = NULL,
                                  sensitivity = NULL, A = NULL,
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
  metric <- match.arg(metric)
  current <- .as_numeric_matrix(current, "current")
  future <- .as_numeric_matrix(future, "future")
  if (!identical(dim(current), dim(future))) {
    stop("current and future must have identical dimensions.", call. = FALSE)
  }
  if (boundary <= 0 || boundary >= 1) {
    stop("boundary must be between 0 and 1.", call. = FALSE)
  }

  occupied_weight <- .reference_weights(occupied, nrow(current),
                                        threshold = occupied_threshold)
  occ <- .positive_reference_indices(occupied_weight)
  scaled <- .standardize_pair(current, future, scale, global_mean, global_sd)
  x0 <- scaled$current
  x1 <- scaled$future

  if (is.null(center)) {
    center <- .extract_slot(cnfa, "mf")
  }
  if (is.null(center)) {
    center <- .weighted_mean(x0, occupied_weight)
  }
  center <- as.numeric(center)
  if (length(center) != ncol(current)) {
    stop("center must have length equal to ncol(current).", call. = FALSE)
  }

  if (is.null(A)) {
    if (is.null(sensitivity)) {
      sensitivity <- .extract_slot(cnfa, "sf")
    }
    if (is.null(sensitivity)) {
      sensitivity <- rep(1, ncol(current))
    }
    A <- niche_metric(sensitivity = sensitivity, cnfa = cnfa, type = metric)
  }
  A <- .validate_niche_metric(A, ncol(current), "A")

  psi0 <- niche_potential(x0, center = center, A = A)
  psi1 <- niche_potential(x1, center = center, A = A)
  radius0 <- niche_radius(psi0)
  radius1 <- niche_radius(psi1)
  amount <- .climate_change_amount(x0, x1, A = A)
  distance_change <- .niche_distance_change(psi0, psi1)
  reconfiguration <- .climate_reconfiguration(amount, distance_change)
  alignment <- .change_alignment(amount, distance_change)
  b_potential <- as.numeric(.weighted_quantile(
    psi0,
    weights = occupied_weight,
    probs = boundary,
    names = FALSE
  ))
  b_radius <- sqrt(pmax(0, b_potential))
  exceed <- boundary_exceedance(psi1, boundary_value = b_potential,
                                scale = "radial")
  perc <- niche_percentile(psi0, psi1, occupied = occupied_weight)
  contrib <- variable_contribution(x0, x1, center = center, A = A)
  mixed <- mixed_variable_response(
    contribution = contrib,
    niche_boundary_exceedance = exceed,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance,
    conflict_ratio = conflict_ratio
  )
  class <- classify_exposure(
    climate_change_amount = amount,
    niche_distance_change = distance_change,
    niche_boundary_exceedance = exceed,
    climate_reconfiguration = reconfiguration,
    contribution = contrib,
    tolerance = tolerance,
    tolerance_quantile = tolerance_quantile,
    stable_climate_change = stable_climate_change,
    stable_quantile = stable_quantile,
    stable_reconfiguration = stable_reconfiguration,
    stable_reconfiguration_quantile = stable_reconfiguration_quantile,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance,
    conflict_ratio = conflict_ratio
  )
  classification_settings <- attr(class, "classification_settings")
  attr(class, "classification_settings") <- NULL
  descriptors <- .exposure_descriptors(
    niche_distance_change = distance_change,
    niche_boundary_exceedance = exceed,
    tolerance = classification_settings$tolerance,
    boundary_exceedance_tolerance =
      classification_settings$boundary_exceedance_tolerance
  )

  out <- list(
    call = match.call(),
    current = x0,
    future = x1,
    occupied = occ,
    occupied_weight = occupied_weight,
    reference_weight = occupied_weight,
    center = center,
    A = A,
    psi_current = psi0,
    psi_future = psi1,
    niche_radius_current = radius0,
    niche_radius_future = radius1,
    climate_change_amount = amount,
    niche_distance_change = distance_change,
    climate_reconfiguration = reconfiguration,
    composition_change = reconfiguration,
    change_alignment = alignment,
    niche_boundary_exceedance = exceed,
    outside_niche_exceedance = exceed,
    boundary_quantile = boundary,
    boundary_value = b_radius,
    boundary_potential = b_potential,
    boundary_radius = b_radius,
    niche_percentile = perc,
    variable_contribution = contrib,
    mixed_variable_response = mixed,
    radial_direction = descriptors$radial_direction,
    boundary_status = descriptors$boundary_status,
    classification = class,
    classification_settings = classification_settings,
    standardization = list(center = scaled$center, scale = scaled$scale)
  )
  class(out) <- "climniche_fit"
  out
}

#' Niche potential
#'
#' @param x Standardised climate matrix.
#' @param center Realised niche centre.
#' @param A Niche metric matrix.
#'
#' @return Numeric vector of quadratic niche displacement values.
#' @export
niche_potential <- function(x, center, A) {
  x <- .as_numeric_matrix(x, "x")
  z <- sweep(x, 2L, center, "-")
  .quad_form_rows(z, A)
}

#' Niche radius
#'
#' @param psi Numeric niche-potential values.
#'
#' @return Numeric vector in sensitivity weighted climate distance units.
#' @export
niche_radius <- function(psi) {
  sqrt(pmax(0, psi))
}

.climate_change_amount <- function(current, future, A) {
  current <- .as_numeric_matrix(current, "current")
  future <- .as_numeric_matrix(future, "future")
  sqrt(pmax(0, .quad_form_rows(future - current, A)))
}

.climate_reconfiguration <- function(climate_change_amount,
                                     niche_distance_change) {
  sqrt(pmax(0, climate_change_amount^2 - niche_distance_change^2))
}

.composition_change <- .climate_reconfiguration

.change_alignment <- function(climate_change_amount, niche_distance_change) {
  out <- rep(NA_real_, length(climate_change_amount))
  ok <- is.finite(climate_change_amount) &
    climate_change_amount > sqrt(.Machine$double.eps)
  out[ok] <- niche_distance_change[ok] / climate_change_amount[ok]
  out
}

.niche_distance_change <- function(psi_current, psi_future) {
  niche_radius(psi_future) - niche_radius(psi_current)
}

#' Niche Boundary Exceedance
#'
#' @param psi_future Future niche potential.
#' @param boundary_value Empirical boundary of the current realised niche in
#'   potential units.
#' @param scale `"radial"` returns exceedance beyond the niche boundary distance;
#'   `"potential"` returns exceedance beyond squared niche potential.
#'
#' @return Numeric vector.
#' @export
boundary_exceedance <- function(psi_future, boundary_value,
                                scale = c("radial", "potential")) {
  scale <- match.arg(scale)
  if (scale == "radial") {
    return(pmax(0, niche_radius(psi_future) - niche_radius(boundary_value)))
  }
  pmax(0, psi_future - boundary_value)
}

#' Niche percentile shift
#'
#' @param psi_current Current niche potential for all cells.
#' @param psi_future Future niche potential for all cells.
#' @param occupied Current reference weights or indices used to define the
#'   reference CDF.
#'
#' @return Data frame with current, future, and delta percentiles.
#' @export
niche_percentile <- function(psi_current, psi_future, occupied) {
  weights <- .reference_weights(occupied, length(psi_current), threshold = NULL)
  current <- .weighted_ecdf_values(psi_current, psi_current, weights)
  future <- .weighted_ecdf_values(psi_current, psi_future, weights)
  data.frame(current = current, future = future, delta = future - current)
}

#' Variable contribution to change in niche potential
#'
#' @param current Current standardized climate matrix.
#' @param future Future standardized climate matrix.
#' @param center Realised niche centre.
#' @param A Niche metric matrix.
#'
#' @return Matrix whose rows sum to the change in niche potential.
#' @export
variable_contribution <- function(current, future, center, A) {
  current <- .as_numeric_matrix(current, "current")
  future <- .as_numeric_matrix(future, "future")
  A <- .validate_niche_metric(A, ncol(current), "A")
  c0 <- sweep(current, 2L, center, "-")
  c1 <- sweep(future, 2L, center, "-")
  out <- c1 * (c1 %*% A) - c0 * (c0 %*% A)
  colnames(out) <- colnames(current)
  out
}
