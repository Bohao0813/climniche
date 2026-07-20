#' Fit niche relative climate exposure
#'
#' @param current Numeric matrix or data frame of current climate values.
#' @param future Numeric matrix or data frame of future climate values, with the
#'   same rows and variables as `current`. Complete names are matched before
#'   fitting.
#' @param occupied NULL, logical vector, row indices, or a numeric vector with
#'   one value per row. Numeric vectors of length `nrow(current)` are treated as
#'   continuous reference weights, including SDM suitability values.
#' @param occupied_threshold Threshold used when `occupied` is a numeric vector
#'   with one value per row.
#' @param cnfa Optional CENFA `cnfa` object. When supplied, its `mf` and `sf`
#'   slots are used unless `center`, `sensitivity`, or `A` are provided.
#' @param center Optional realised niche centre in the fitted climate space.
#' @param sensitivity Optional climate-variable sensitivity weights.
#' @param A Optional niche metric matrix.
#' @param metric `"diag"` or `"factor"`. Used only when `A` is missing.
#' @param boundary Quantile for the empirical occupied niche boundary.
#' @param scale Logical. If TRUE, center and scale current and future values
#'   using means and standard deviations of current values.
#' @param preprocess Logical. If TRUE, remove near-zero variance variables and
#'   highly correlated variables before scaling and metric fitting.
#' @param preprocess_correlation Maximum absolute current-climate correlation
#'   retained during preprocessing.
#' @param preprocess_min_sd Minimum current-climate standard deviation retained
#'   during preprocessing.
#' @param global_mean Optional means used for centering.
#' @param global_sd Optional standard deviations used for scaling.
#' @param tolerance Optional Niche Distance Shift tolerance.
#' @param tolerance_quantile Quantile of absolute Niche Distance Shift used
#'   when `tolerance = NULL`.
#' @param boundary_exceedance_tolerance Tolerance for deciding whether future
#'   climate exceeds the empirical niche boundary.
#'
#' @return An object of class `climniche_fit`.
#' @noRd
.fit_climniche_matrix <- function(current, future, occupied = NULL,
                                  occupied_threshold = NULL,
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
  metric <- match.arg(metric)
  metric_type <- if (is.null(A)) metric else "user"
  scale <- .check_flag(scale, "scale")
  preprocess <- .check_flag(preprocess, "preprocess")
  boundary <- .check_open_probability(boundary, "boundary")
  current <- .as_numeric_matrix(current, "current")
  future <- .as_numeric_matrix(future, "future")
  aligned <- .align_climate_pair(current, future)
  current <- aligned$current
  future <- aligned$future

  occupied_weight <- .reference_weights(
    occupied,
    nrow(current),
    threshold = occupied_threshold,
    row_names = rownames(current)
  )
  occ <- .positive_reference_indices(occupied_weight)
  p_original <- ncol(current)
  preprocessed <- .preprocess_climate_pair(
    current = current,
    future = future,
    preprocess = preprocess,
    correlation = preprocess_correlation,
    min_sd = preprocess_min_sd
  )
  current <- preprocessed$current
  future <- preprocessed$future
  keep <- preprocessed$keep
  preprocessing_record <- preprocessed[
    c("keep", "original_variables", "retained_variables",
      "removed_variables", "settings")
  ]
  variables <- preprocessed$original_variables
  center <- .subset_fit_vector(center, keep, p_original, "center", variables)
  sensitivity <- .subset_fit_vector(sensitivity, keep, p_original,
                                    "sensitivity", variables)
  global_mean <- .subset_fit_vector(global_mean, keep, p_original,
                                    "global_mean", variables)
  global_sd <- .subset_fit_vector(global_sd, keep, p_original, "global_sd",
                                  variables)
  A <- .subset_fit_matrix(A, keep, p_original, "A", variables)
  cnfa <- .subset_cnfa_object(cnfa, keep, p_original, variables)
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
  if (length(center) != ncol(current) || any(!is.finite(center))) {
    stop("center must contain one finite value per retained climate variable.",
         call. = FALSE)
  }
  names(center) <- colnames(current)

  if (is.null(A)) {
    if (identical(metric, "diag")) {
      if (is.null(sensitivity)) {
        sensitivity <- .extract_slot(cnfa, "sf")
      }
      if (is.null(sensitivity)) {
        sensitivity <- rep(1, ncol(current))
      }
      if (length(sensitivity) != ncol(current)) {
        stop("sensitivity must contain one value per retained climate variable.",
             call. = FALSE)
      }
    }
    A <- niche_metric(sensitivity = sensitivity, cnfa = cnfa, type = metric)
  }
  A <- .validate_niche_metric(A, ncol(current), "A")
  dimnames(A) <- list(colnames(current), colnames(current))
  sensitivity_weights <- if (identical(metric_type, "diag")) {
    sensitivity <- as.numeric(sensitivity)
    sensitivity / mean(sensitivity)
  } else {
    NULL
  }
  if (!is.null(sensitivity_weights)) {
    names(sensitivity_weights) <- colnames(current)
  }

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
  descriptors <- .exposure_descriptors(
    niche_distance_change = distance_change,
    niche_boundary_exceedance = exceed,
    tolerance = tolerance,
    tolerance_quantile = tolerance_quantile,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance
  )
  descriptor_settings <- attr(descriptors, "descriptor_settings")
  attr(descriptors, "descriptor_settings") <- NULL

  out <- list(
    call = match.call(),
    current = x0,
    future = x1,
    occupied = occ,
    occupied_weight = occupied_weight,
    reference_weight = occupied_weight,
    center = center,
    A = A,
    metric_type = metric_type,
    sensitivity_weights = sensitivity_weights,
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
    boundary_distance = b_radius,
    boundary_potential = b_potential,
    boundary_radius = b_radius,
    niche_percentile = perc,
    variable_contribution = contrib,
    radial_direction = descriptors$radial_direction,
    boundary_status = descriptors$boundary_status,
    descriptor_settings = descriptor_settings,
    threshold_settings = descriptor_settings,
    preprocessing = preprocessing_record,
    standardization = list(
      enabled = scale,
      mean = scaled$center,
      sd = scaled$scale,
      center = scaled$center,
      scale = scaled$scale
    )
  )
  class(out) <- "climniche_fit"
  out
}

#' Niche potential
#'
#' @param x Climate matrix in the fitted climate space.
#' @param center Realised niche centre.
#' @param A Niche metric matrix.
#'
#' @return Numeric vector of quadratic niche displacement values.
#' @export
niche_potential <- function(x, center, A) {
  x <- .as_numeric_matrix(x, "x")
  aligned <- .align_metric_inputs(x, center, A)
  center <- aligned$center
  A <- aligned$A
  if (length(center) != ncol(x) || any(!is.finite(center))) {
    stop("center must contain one finite value per column of x.",
         call. = FALSE)
  }
  A <- .validate_niche_metric(A, ncol(x), "A")
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
#'
#' @examples
#' psi <- c(0.5, 1, 2)
#' boundary_exceedance(psi, boundary_value = 1)
#' @export
boundary_exceedance <- function(psi_future, boundary_value,
                                scale = c("radial", "potential")) {
  scale <- match.arg(scale)
  psi_future <- as.numeric(psi_future)
  boundary_value <- .check_finite_scalar(boundary_value, "boundary_value")
  if (any(!is.finite(psi_future))) {
    stop("psi_future must contain finite values.", call. = FALSE)
  }
  if (boundary_value < 0) {
    stop("boundary_value must be non-negative.", call. = FALSE)
  }
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
  if (length(psi_current) != length(psi_future)) {
    stop("psi_current and psi_future must have the same length.",
         call. = FALSE)
  }
  weights <- .reference_weights(occupied, length(psi_current), threshold = NULL)
  current <- .weighted_ecdf_values(psi_current, psi_current, weights)
  future <- .weighted_ecdf_values(psi_current, psi_future, weights)
  data.frame(current = current, future = future, delta = future - current)
}

#' Variable contribution to change in niche potential
#'
#' @param current Current climate matrix in the fitted climate space.
#' @param future Future climate matrix in the fitted climate space.
#' @param center Realised niche centre.
#' @param A Niche metric matrix.
#'
#' @return Matrix whose rows sum to the change in niche potential.
#' @export
variable_contribution <- function(current, future, center, A) {
  current <- .as_numeric_matrix(current, "current")
  future <- .as_numeric_matrix(future, "future")
  aligned <- .align_climate_pair(current, future)
  current <- aligned$current
  future <- aligned$future
  metric_inputs <- .align_metric_inputs(current, center, A)
  center <- metric_inputs$center
  A <- metric_inputs$A
  if (length(center) != ncol(current) || any(!is.finite(center))) {
    stop("center must contain one finite value per climate variable.",
         call. = FALSE)
  }
  A <- .validate_niche_metric(A, ncol(current), "A")
  c0 <- sweep(current, 2L, center, "-")
  c1 <- sweep(future, 2L, center, "-")
  out <- c1 * (c1 %*% A) - c0 * (c0 %*% A)
  colnames(out) <- colnames(current)
  out
}
