#' Fit a current climatic niche reference
#'
#' Estimates the preprocessing, standardisation, realised niche centre,
#' climatic weighting matrix and empirical niche boundary once. The resulting
#' object can be reused with [project_climniche()] for several future periods or
#' climate models.
#'
#' @param current Numeric matrix or data frame of current climatic conditions.
#'   Rows are cells or sites and columns are climatic variables.
#' @param occupied Reference cells or weights, following the rules used by
#'   [fit_climniche()].
#' @param occupied_threshold Optional cutoff for numeric reference weights.
#' @param cnfa Optional compatible CENFA object.
#' @param center Optional realised niche centre in the fitted climatic space.
#' @param sensitivity Optional non-negative variable sensitivity weights.
#' @param A Optional climatic weighting matrix.
#' @param metric Method used to build `A` when it is not supplied.
#' @param boundary Weighted quantile defining the empirical niche boundary.
#' @param scale If `TRUE`, standardise retained variables using current-climate
#'   means and standard deviations.
#' @param preprocess If `TRUE`, remove near-zero variance and highly correlated
#'   variables before fitting.
#' @param preprocess_correlation Maximum absolute correlation retained during
#'   preprocessing.
#' @param preprocess_min_sd Minimum current-climate standard deviation retained
#'   during preprocessing.
#' @param global_mean,global_sd Optional centring and scaling values.
#'
#' @return A `climniche_reference` object.
#'
#' @details
#' The reference object fixes the climatic space used for all subsequent
#' projections. Future conditions do not alter the centre, weighting matrix,
#' standardisation or empirical boundary.
#'
#' @examples
#' sim <- simulate_climniche(n = 200, p = 6, seed = 4)
#' reference <- fit_climniche_reference(
#'   sim$current,
#'   occupied = sim$occupied,
#'   sensitivity = sim$sensitivity
#' )
#' reference
#' @export
fit_climniche_reference <- function(
    current,
    occupied = NULL,
    occupied_threshold = NULL,
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
    global_sd = NULL) {
  metric <- match.arg(metric)
  fit <- .fit_climniche_matrix(
    current = current,
    future = current,
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
    tolerance = 0,
    tolerance_quantile = 0.10,
    boundary_exceedance_tolerance = 0
  )
  out <- .reference_from_fit(fit)
  out$call <- match.call()
  out
}

.reference_from_fit <- function(fit) {
  if (!inherits(fit, "climniche_fit")) {
    stop("fit must be a climniche_fit object.", call. = FALSE)
  }
  reference_weight <- .fit_reference_weights(fit)
  row_names <- rownames(fit$current)
  if (.names_are_complete(row_names)) {
    names(reference_weight) <- row_names
  }
  out <- list(
    call = fit$call,
    current = fit$current,
    reference_weight = reference_weight,
    occupied_weight = reference_weight,
    occupied = .positive_reference_indices(reference_weight),
    center = fit$center,
    A = fit$A,
    metric_type = fit$metric_type,
    sensitivity_weights = fit$sensitivity_weights,
    psi_reference = fit$psi_current,
    niche_radius_reference = fit$niche_radius_current,
    boundary_quantile = fit$boundary_quantile,
    boundary_value = fit$boundary_radius,
    boundary_distance = fit$boundary_radius,
    boundary_potential = fit$boundary_potential,
    boundary_radius = fit$boundary_radius,
    preprocessing = fit$preprocessing,
    standardization = fit$standardization
  )
  class(out) <- "climniche_reference"
  out
}

#' Project climatic conditions onto a fitted niche reference
#'
#' @param reference A `climniche_reference` object.
#' @param future Numeric matrix or data frame of projected climatic conditions.
#' @param current Optional current climatic conditions for the projected rows.
#'   If omitted, the current rows stored in `reference` are used.
#' @param occupied Optional reference weights for summaries of the projected
#'   rows. When `current` is omitted, the fitted reference weights are reused.
#'   With a supplied `current`, named rows are matched to the reference;
#'   otherwise projected rows receive equal weight.
#' @param occupied_threshold Optional cutoff for numeric `occupied` weights.
#' @param tolerance Optional tolerance around zero for Niche Distance Shift.
#' @param tolerance_quantile Quantile of absolute Niche Distance Shift used
#'   when `tolerance = NULL`.
#' @param boundary_exceedance_tolerance Non-negative tolerance used by the
#'   boundary status descriptor.
#'
#' @return A `climniche_fit` object using the fixed reference.
#'
#' @details
#' Let \eqn{A} be the fitted weighting matrix and let \eqn{\mu} be the realised
#' niche centre. In the transformed climatic space, write the current and
#' future centred vectors as \eqn{z_0} and \eqn{z_1}, with lengths \eqn{r_0}
#' and \eqn{r_1} and angle \eqn{\theta}. Climatic Reconfiguration satisfies
#' \deqn{C_i^2 = 2 r_{0i} r_{1i} (1 - \cos(\theta_i)).}
#' It therefore combines angular change with current and future niche
#' distances and is calculated rather than fitted independently.
#'
#' @examples
#' sim <- simulate_climniche(n = 200, p = 6, seed = 4)
#' reference <- fit_climniche_reference(
#'   sim$current,
#'   occupied = sim$occupied,
#'   sensitivity = sim$sensitivity
#' )
#' fit <- project_climniche(reference, sim$future_away)
#' climniche_summary(fit)
#' @export
project_climniche <- function(reference, future, current = NULL,
                              occupied = NULL,
                              occupied_threshold = NULL,
                              tolerance = NULL,
                              tolerance_quantile = 0.10,
                              boundary_exceedance_tolerance = 0) {
  if (!inherits(reference, "climniche_reference")) {
    stop("reference must be a climniche_reference object.", call. = FALSE)
  }

  reuse_reference <- is.null(current)
  if (reuse_reference) {
    x0 <- reference$current
    x1 <- .transform_with_reference(future, reference, "future")
    aligned <- .align_climate_pair(x0, x1)
  } else {
    current <- .as_numeric_matrix(current, "current")
    future <- .as_numeric_matrix(future, "future")
    raw <- .align_climate_pair(current, future)
    x0 <- .transform_with_reference(raw$current, reference, "current")
    x1 <- .transform_with_reference(raw$future, reference, "future")
    aligned <- .align_climate_pair(x0, x1)
  }
  x0 <- aligned$current
  x1 <- aligned$future

  occupied_weight <- .projection_weights(
    occupied = occupied,
    occupied_threshold = occupied_threshold,
    reference = reference,
    row_names = rownames(x0),
    n = nrow(x0),
    reuse_reference = reuse_reference
  )
  occ <- .positive_reference_indices(occupied_weight)

  psi0 <- niche_potential(x0, center = reference$center, A = reference$A)
  psi1 <- niche_potential(x1, center = reference$center, A = reference$A)
  radius0 <- niche_radius(psi0)
  radius1 <- niche_radius(psi1)
  amount <- .climate_change_amount(x0, x1, A = reference$A)
  distance_change <- .niche_distance_change(psi0, psi1)
  reconfiguration <- .climate_reconfiguration(amount, distance_change)
  alignment <- .change_alignment(amount, distance_change)
  exceed <- boundary_exceedance(
    psi1,
    boundary_value = reference$boundary_potential,
    scale = "radial"
  )
  percentile_current <- .weighted_ecdf_values(
    reference$psi_reference,
    psi0,
    reference$reference_weight
  )
  percentile_future <- .weighted_ecdf_values(
    reference$psi_reference,
    psi1,
    reference$reference_weight
  )
  perc <- data.frame(
    current = percentile_current,
    future = percentile_future,
    delta = percentile_future - percentile_current
  )
  contrib <- variable_contribution(
    x0,
    x1,
    center = reference$center,
    A = reference$A
  )
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
    center = reference$center,
    A = reference$A,
    metric_type = reference$metric_type,
    sensitivity_weights = reference$sensitivity_weights,
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
    boundary_quantile = reference$boundary_quantile,
    boundary_value = reference$boundary_radius,
    boundary_distance = reference$boundary_radius,
    boundary_potential = reference$boundary_potential,
    boundary_radius = reference$boundary_radius,
    niche_percentile = perc,
    variable_contribution = contrib,
    radial_direction = descriptors$radial_direction,
    boundary_status = descriptors$boundary_status,
    descriptor_settings = descriptor_settings,
    threshold_settings = descriptor_settings,
    preprocessing = reference$preprocessing,
    standardization = reference$standardization,
    reference_call = reference$call
  )
  class(out) <- "climniche_fit"
  out
}

.transform_with_reference <- function(x, reference, arg) {
  x <- .as_numeric_matrix(x, arg)
  original <- reference$preprocessing$original_variables
  retained <- reference$preprocessing$retained_variables
  variables <- colnames(x)

  if (.names_are_partial(variables)) {
    stop(arg, " variable names must be complete or omitted.", call. = FALSE)
  }
  if (.names_are_complete(variables)) {
    .check_unique_names(variables, paste(arg, "variable names"))
    target <- if (setequal(variables, original)) {
      original
    } else if (setequal(variables, retained)) {
      retained
    } else {
      stop(arg, " variables must match the fitted reference.", call. = FALSE)
    }
    x <- x[, match(target, variables), drop = FALSE]
  } else if (ncol(x) == length(original)) {
    colnames(x) <- original
  } else if (ncol(x) == length(retained)) {
    colnames(x) <- retained
  } else {
    stop(arg, " must contain the original or retained reference variables.",
         call. = FALSE)
  }

  if (identical(colnames(x), original)) {
    x <- x[, retained, drop = FALSE]
  }
  x <- x[, retained, drop = FALSE]
  if (any(!is.finite(x))) {
    stop(arg, " must contain finite values.", call. = FALSE)
  }

  centre <- reference$standardization$center
  scale_value <- reference$standardization$scale
  centre <- centre[match(retained, names(centre))]
  scale_value <- scale_value[match(retained, names(scale_value))]
  if (any(!is.finite(centre)) || any(!is.finite(scale_value)) ||
      any(scale_value <= 0)) {
    stop("reference contains invalid standardisation values.", call. = FALSE)
  }
  sweep(sweep(x, 2L, centre, "-"), 2L, scale_value, "/")
}

.projection_weights <- function(occupied, occupied_threshold, reference,
                                row_names, n, reuse_reference = FALSE) {
  if (!is.null(occupied)) {
    return(.reference_weights(
      occupied,
      n,
      threshold = occupied_threshold,
      row_names = row_names
    ))
  }
  weights <- reference$reference_weight
  if (length(weights) != n && reuse_reference) {
    return(rep(1, n))
  }
  weight_names <- names(weights)
  if (.names_are_complete(weight_names) && .names_are_complete(row_names)) {
    if (!setequal(weight_names, row_names)) {
      return(rep(1, n))
    }
    return(unname(weights[match(row_names, weight_names)]))
  }
  if (reuse_reference && length(weights) == n) {
    return(unname(weights))
  }
  rep(1, n)
}

#' @export
print.climniche_reference <- function(x, ...) {
  cat("Current climatic niche reference\n")
  cat("Reference rows:", length(x$reference_weight), "\n")
  cat("Positive reference weights:", sum(x$reference_weight > 0), "\n")
  cat("Retained variables:", ncol(x$current), "\n")
  cat("Boundary quantile:", format(x$boundary_quantile), "\n")
  cat("Boundary distance:", format(x$boundary_distance), "\n")
  invisible(x)
}
