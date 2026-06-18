#' @export
summary.climniche_fit <- function(object, ...) {
  descriptors <- .fit_exposure_descriptors(object)
  out <- list(
    n = nrow(object$current),
    p = ncol(object$current),
    occupied_n = length(object$occupied),
    boundary_quantile = object$boundary_quantile,
    boundary_value = object$boundary_value,
    climate_change_amount = base::summary(object$climate_change_amount),
    niche_distance_change = base::summary(object$niche_distance_change),
    climate_reconfiguration =
      base::summary(.fit_metric(object, "climate_reconfiguration")),
    composition_change = base::summary(.fit_metric(object, "composition_change")),
    change_alignment = base::summary(object$change_alignment),
    niche_boundary_exceedance =
      base::summary(.fit_metric(object, "niche_boundary_exceedance")),
    outside_niche_exceedance =
      base::summary(.fit_metric(object, "outside_niche_exceedance")),
    classification_settings = object$classification_settings,
    radial_direction = table(descriptors$radial_direction),
    boundary_status = table(descriptors$boundary_status),
    classification = table(.normalise_class(object$classification))
  )
  class(out) <- "summary.climniche_fit"
  out
}

#' @export
plot.climniche_fit <- function(x, type = c("distance", "boundary", "amount",
                                           "reconfiguration",
                                           "classification"), ...) {
  type <- match.arg(type)
  if (type == "distance") {
    graphics::hist(x$niche_distance_change, main = "Niche Distance Shift",
                   xlab = "Niche Distance Shift", ...)
  } else if (type == "boundary") {
    graphics::hist(.fit_metric(x, "niche_boundary_exceedance"),
                   main = "Niche Boundary Exceedance",
                   xlab = "Niche Boundary Exceedance", ...)
  } else if (type == "amount") {
    graphics::hist(x$climate_change_amount, main = "Climatic Displacement",
                   xlab = "Climatic Displacement", ...)
  } else if (type == "reconfiguration") {
    graphics::hist(.fit_metric(x, "climate_reconfiguration"),
                   main = "Climatic Reconfiguration",
                   xlab = "Climatic Reconfiguration", ...)
  } else {
    graphics::barplot(table(x$classification), las = 2,
                      main = "Derived exposure classes", ...)
  }
  invisible(x)
}
