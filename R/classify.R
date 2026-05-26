#' Classify projected climate change relative to the current niche
#'
#' @param climate_change_amount Numeric climate change amount.
#' @param niche_distance_change Numeric change in distance to the current niche.
#' @param outside_niche_exceedance Numeric exceedance beyond the niche boundary.
#' @param composition_change Optional composition change.
#' @param contribution Optional variable contribution matrix.
#' @param tolerance Optional niche distance change tolerance. If NULL, the 10th
#'   percentile of absolute niche distance change is used.
#' @param stable_quantile Quantile of climate change amount used to mark weak
#'   climate change.
#' @param conflict_ratio Retained for compatibility with the fitting workflow.
#'
#' @return Factor of climniche classes.
classify_exposure <- function(climate_change_amount, niche_distance_change,
                              outside_niche_exceedance,
                              composition_change = NULL,
                              contribution = NULL,
                              tolerance = NULL, stable_quantile = 0.25,
                              conflict_ratio = NULL) {
  if (is.null(tolerance)) {
    tolerance <- as.numeric(stats::quantile(abs(niche_distance_change),
                                            probs = 0.10, names = FALSE,
                                            na.rm = TRUE, type = 8))
  }
  stable_cut <- as.numeric(stats::quantile(
    climate_change_amount,
    probs = stable_quantile,
    names = FALSE,
    na.rm = TRUE,
    type = 8
  ))
  cls <- rep("changed composition, similar distance",
             length(climate_change_amount))
  cls[niche_distance_change < -tolerance] <- "closer to current niche"
  cls[niche_distance_change > tolerance] <- "farther from current niche"
  cls[outside_niche_exceedance > 0] <- "outside current niche boundary"
  cls[climate_change_amount <= stable_cut & outside_niche_exceedance == 0 &
        abs(niche_distance_change) <= tolerance] <- "little climate niche change"

  if (!is.null(composition_change)) {
    low_composition <- composition_change <= stats::quantile(
      composition_change,
      probs = stable_quantile,
      names = FALSE,
      na.rm = TRUE,
      type = 8
    )
    cls[low_composition & climate_change_amount <= stable_cut &
          outside_niche_exceedance == 0 &
          abs(niche_distance_change) <= tolerance] <- "little climate niche change"
  }

  factor(cls, levels = c(
    "little climate niche change",
    "closer to current niche",
    "farther from current niche",
    "outside current niche boundary",
    "changed composition, similar distance"
  ))
}

mixed_variable_response <- function(contribution, outside_niche_exceedance = NULL,
                                    conflict_ratio = 0.25) {
  if (is.null(contribution) || is.null(conflict_ratio)) {
    n <- if (is.null(outside_niche_exceedance)) 0 else length(outside_niche_exceedance)
    return(rep(FALSE, n))
  }
  pos <- rowSums(pmax(contribution, 0))
  neg <- rowSums(abs(pmin(contribution, 0)))
  total <- pos + neg
  mixed <- total > 0 & pmin(pos, neg) / total >= conflict_ratio
  if (!is.null(outside_niche_exceedance)) {
    mixed <- mixed & outside_niche_exceedance == 0
  }
  mixed
}
