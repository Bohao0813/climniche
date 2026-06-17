#' Classify projected climate change relative to the current niche
#'
#' @param climate_change_amount Climatic Displacement: total sensitivity
#'   weighted climatic displacement.
#' @param niche_distance_change Niche Distance Shift: signed change in distance
#'   from the current realised climatic niche centre.
#' @param niche_boundary_exceedance Niche Boundary Exceedance: positive excess
#'   of future niche distance beyond the empirical boundary of the current
#'   realised climatic niche.
#' @param climate_reconfiguration Climatic Reconfiguration: non radial
#'   component of climatic displacement not captured by change in distance to
#'   the niche centre.
#' @param contribution Optional variable contribution matrix.
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
#' @param conflict_ratio Retained for compatibility with the fitting workflow.
#' @param outside_niche_exceedance,composition_change Backwards-compatible
#'   aliases.
#'
#' @return Factor of climniche classes.
#' @noRd
classify_exposure <- function(climate_change_amount, niche_distance_change,
                              niche_boundary_exceedance = NULL,
                              climate_reconfiguration = NULL,
                              outside_niche_exceedance = NULL,
                              composition_change = NULL,
                              contribution = NULL,
                              tolerance = NULL,
                              tolerance_quantile = 0.10,
                              stable_climate_change = NULL,
                              stable_quantile = 0.25,
                              stable_reconfiguration = NULL,
                              stable_reconfiguration_quantile = 0.25,
                              boundary_exceedance_tolerance = 0,
                              conflict_ratio = NULL) {
  if (is.null(niche_boundary_exceedance)) {
    niche_boundary_exceedance <- outside_niche_exceedance
  }
  if (is.null(climate_reconfiguration)) {
    climate_reconfiguration <- composition_change
  }
  if (is.null(niche_boundary_exceedance)) {
    stop("niche_boundary_exceedance must be supplied.", call. = FALSE)
  }
  if (is.null(climate_reconfiguration)) {
    climate_reconfiguration <- .climate_reconfiguration(
      climate_change_amount, niche_distance_change
    )
  }
  .check_quantile <- function(x, name) {
    x <- as.numeric(x)[1]
    if (!is.finite(x) || x < 0 || x > 1) {
      stop(name, " must be a finite value between 0 and 1.", call. = FALSE)
    }
    x
  }
  tolerance_quantile <- .check_quantile(tolerance_quantile,
                                        "tolerance_quantile")
  stable_quantile <- .check_quantile(stable_quantile, "stable_quantile")
  stable_reconfiguration_quantile <- .check_quantile(
    stable_reconfiguration_quantile,
    "stable_reconfiguration_quantile"
  )
  boundary_exceedance_tolerance <- as.numeric(boundary_exceedance_tolerance)[1]
  if (!is.finite(boundary_exceedance_tolerance) ||
      boundary_exceedance_tolerance < 0) {
    stop("boundary_exceedance_tolerance must be a finite non-negative number.",
         call. = FALSE)
  }
  if (is.null(tolerance)) {
    tolerance <- as.numeric(stats::quantile(abs(niche_distance_change),
                                            probs = tolerance_quantile,
                                            names = FALSE,
                                            na.rm = TRUE, type = 8))
  } else {
    tolerance <- as.numeric(tolerance)[1]
  }
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a finite non-negative number.", call. = FALSE)
  }
  if (is.null(stable_climate_change)) {
    stable_climate_change <- as.numeric(stats::quantile(
      climate_change_amount,
      probs = stable_quantile,
      names = FALSE,
      na.rm = TRUE,
      type = 8
    ))
  } else {
    stable_climate_change <- as.numeric(stable_climate_change)[1]
  }
  if (!is.finite(stable_climate_change) || stable_climate_change < 0) {
    stop("stable_climate_change must be a finite non-negative number.",
         call. = FALSE)
  }
  if (is.null(stable_reconfiguration)) {
    stable_reconfiguration <- as.numeric(stats::quantile(
      climate_reconfiguration,
      probs = stable_reconfiguration_quantile,
      names = FALSE,
      na.rm = TRUE,
      type = 8
    ))
  } else {
    stable_reconfiguration <- as.numeric(stable_reconfiguration)[1]
  }
  if (!is.finite(stable_reconfiguration) || stable_reconfiguration < 0) {
    stop("stable_reconfiguration must be a finite non-negative number.",
         call. = FALSE)
  }

  cls <- rep("Climatic Reconfiguration with limited Niche Distance Shift",
             length(climate_change_amount))
  cls[niche_distance_change < -tolerance] <- "Closer to current niche"
  cls[niche_distance_change > tolerance] <- "Farther from current niche"
  exceeds <- niche_boundary_exceedance > boundary_exceedance_tolerance
  cls[exceeds] <- "Outside current niche boundary"
  stable <- climate_change_amount <= stable_climate_change &
    climate_reconfiguration <= stable_reconfiguration &
    !exceeds &
    abs(niche_distance_change) <= tolerance
  cls[stable] <- "Limited niche relative change"

  out <- factor(cls, levels = .class_level_names())
  attr(out, "classification_settings") <- list(
    tolerance = tolerance,
    tolerance_quantile = tolerance_quantile,
    stable_climate_change = stable_climate_change,
    stable_quantile = stable_quantile,
    stable_reconfiguration = stable_reconfiguration,
    stable_reconfiguration_quantile = stable_reconfiguration_quantile,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance,
    conflict_ratio = conflict_ratio
  )
  out
}

mixed_variable_response <- function(contribution, niche_boundary_exceedance = NULL,
                                    outside_niche_exceedance = NULL,
                                    boundary_exceedance_tolerance = 0,
                                    conflict_ratio = 0.25) {
  if (is.null(niche_boundary_exceedance)) {
    niche_boundary_exceedance <- outside_niche_exceedance
  }
  if (is.null(contribution) || is.null(conflict_ratio)) {
    n <- if (is.null(niche_boundary_exceedance)) {
      0
    } else {
      length(niche_boundary_exceedance)
    }
    return(rep(FALSE, n))
  }
  pos <- rowSums(pmax(contribution, 0))
  neg <- rowSums(abs(pmin(contribution, 0)))
  total <- pos + neg
  mixed <- total > 0 & pmin(pos, neg) / total >= conflict_ratio
  if (!is.null(niche_boundary_exceedance)) {
    mixed <- mixed & niche_boundary_exceedance <= boundary_exceedance_tolerance
  }
  mixed
}
