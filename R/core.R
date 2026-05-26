#' Fit niche climate exposure
#'
#' @param current Numeric matrix or data frame of current climate values.
#' @param future Numeric matrix or data frame of future climate values, with the
#'   same dimensions and variables as `current`.
#' @param occupied NULL, logical vector, row indices, or a numeric vector with
#'   one value per row. Numeric vectors of length `nrow(current)` are treated as
#'   binary or continuous occurrence, range, or SDM suitability values.
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
#' @param conflict_ratio Minimum minority-sign variable-contribution share used
#'   to mark mixed variable responses. Set to NULL to disable this flag.
#'
#' @return An object of class `climniche_fit`.
.fit_climniche_matrix <- function(current, future, occupied = NULL,
                                  occupied_threshold = 0,
                                  cnfa = NULL, center = NULL,
                                  sensitivity = NULL, A = NULL,
                                  metric = c("diag", "factor"),
                                  boundary = 0.95, scale = TRUE,
                                  global_mean = NULL, global_sd = NULL,
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

  occ <- .occupied_index(occupied, nrow(current),
                         threshold = occupied_threshold)
  scaled <- .standardize_pair(current, future, scale, global_mean, global_sd)
  x0 <- scaled$current
  x1 <- scaled$future

  if (is.null(center)) {
    center <- .extract_slot(cnfa, "mf")
  }
  if (is.null(center)) {
    center <- colMeans(x0[occ, , drop = FALSE])
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
  composition <- .composition_change(amount, distance_change)
  alignment <- .change_alignment(amount, distance_change)
  b_potential <- as.numeric(stats::quantile(psi0[occ], probs = boundary,
                                            names = FALSE, na.rm = TRUE,
                                            type = 8))
  b_radius <- sqrt(pmax(0, b_potential))
  exceed <- boundary_exceedance(psi1, boundary_value = b_potential,
                                scale = "radial")
  perc <- niche_percentile(psi0, psi1, occupied = occ)
  contrib <- variable_contribution(x0, x1, center = center, A = A)
  mixed <- mixed_variable_response(
    contribution = contrib,
    outside_niche_exceedance = exceed,
    conflict_ratio = conflict_ratio
  )
  class <- classify_exposure(
    climate_change_amount = amount,
    niche_distance_change = distance_change,
    outside_niche_exceedance = exceed,
    composition_change = composition,
    contribution = contrib,
    conflict_ratio = conflict_ratio
  )

  out <- list(
    call = match.call(),
    current = x0,
    future = x1,
    occupied = occ,
    center = center,
    A = A,
    psi_current = psi0,
    psi_future = psi1,
    niche_radius_current = radius0,
    niche_radius_future = radius1,
    climate_change_amount = amount,
    niche_distance_change = distance_change,
    composition_change = composition,
    change_alignment = alignment,
    outside_niche_exceedance = exceed,
    boundary_quantile = boundary,
    boundary_value = b_radius,
    boundary_potential = b_potential,
    boundary_radius = b_radius,
    niche_percentile = perc,
    variable_contribution = contrib,
    mixed_variable_response = mixed,
    classification = class,
    standardization = list(center = scaled$center, scale = scaled$scale)
  )
  class(out) <- "climniche_fit"
  out
}

#' Niche potential
#'
#' @param x Standardized climate matrix.
#' @param center Realized niche center.
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
#' @return Numeric vector in sensitivity-weighted climate-distance units.
#' @export
niche_radius <- function(psi) {
  sqrt(pmax(0, psi))
}

.climate_change_amount <- function(current, future, A) {
  current <- .as_numeric_matrix(current, "current")
  future <- .as_numeric_matrix(future, "future")
  sqrt(pmax(0, .quad_form_rows(future - current, A)))
}

.composition_change <- function(climate_change_amount, niche_distance_change) {
  sqrt(pmax(0, climate_change_amount^2 - niche_distance_change^2))
}

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

#' Niche boundary exceedance
#'
#' @param psi_future Future niche potential.
#' @param boundary_value Empirical occupied-niche boundary in potential units.
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
#' @param occupied Current occurrence indices used to define the reference CDF.
#'
#' @return Data frame with current, future, and delta percentiles.
#' @export
niche_percentile <- function(psi_current, psi_future, occupied) {
  f <- stats::ecdf(psi_current[occupied])
  current <- f(psi_current)
  future <- f(psi_future)
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
