.radial_direction_levels <- function() {
  c(
    "Toward realised niche centre",
    "Limited Niche Distance Shift",
    "Away from realised niche centre"
  )
}

.boundary_status_levels <- function() {
  c(
    "Within empirical niche boundary",
    "Beyond empirical niche boundary"
  )
}

.exposure_descriptors <- function(niche_distance_change,
                                  niche_boundary_exceedance,
                                  tolerance = NULL,
                                  tolerance_quantile = 0.10,
                                  boundary_exceedance_tolerance = 0) {
  tolerance_quantile <- .check_probability(tolerance_quantile,
                                           "tolerance_quantile")
  if (is.null(tolerance)) {
    tolerance <- as.numeric(stats::quantile(
      abs(niche_distance_change),
      probs = tolerance_quantile,
      names = FALSE,
      na.rm = TRUE,
      type = 8
    ))
  } else {
    tolerance <- as.numeric(tolerance)[1]
  }
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a finite non-negative number.", call. = FALSE)
  }

  boundary_exceedance_tolerance <- as.numeric(boundary_exceedance_tolerance)[1]
  if (!is.finite(boundary_exceedance_tolerance) ||
      boundary_exceedance_tolerance < 0) {
    stop("boundary_exceedance_tolerance must be a finite non-negative number.",
         call. = FALSE)
  }

  radial_direction <- rep("Limited Niche Distance Shift",
                          length(niche_distance_change))
  radial_direction[niche_distance_change < -tolerance] <-
    "Toward realised niche centre"
  radial_direction[niche_distance_change > tolerance] <-
    "Away from realised niche centre"

  boundary_status <- ifelse(
    niche_boundary_exceedance > boundary_exceedance_tolerance,
    "Beyond empirical niche boundary",
    "Within empirical niche boundary"
  )

  out <- list(
    radial_direction = factor(
      radial_direction,
      levels = .radial_direction_levels()
    ),
    boundary_status = factor(
      boundary_status,
      levels = .boundary_status_levels()
    )
  )
  attr(out, "descriptor_settings") <- list(
    tolerance = tolerance,
    tolerance_quantile = tolerance_quantile,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance
  )
  out
}
