#' @export
summary.climniche_fit <- function(object, ...) {
  descriptors <- .fit_exposure_descriptors(object)
  boundary_distance <- object$boundary_distance %||%
    object$boundary_radius %||% object$boundary_value
  out <- list(
    n = nrow(object$current),
    p = ncol(object$current),
    occupied_n = length(object$occupied),
    reference_weight_sum = sum(.fit_reference_weights(object)),
    metric_type = object$metric_type %||% "unspecified",
    boundary_quantile = object$boundary_quantile,
    boundary_distance = boundary_distance,
    boundary_value = boundary_distance,
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
    descriptor_settings = object$descriptor_settings,
    radial_direction = table(descriptors$radial_direction),
    boundary_status = table(descriptors$boundary_status)
  )
  class(out) <- "summary.climniche_fit"
  out
}

#' @export
print.summary.climniche_fit <- function(x, digits = max(
  3L, getOption("digits") - 3L
), ...) {
  cat("climniche fit summary\n\n")
  settings <- data.frame(
    Cells = x$n,
    Variables = x$p,
    `Reference cells` = x$occupied_n,
    `Reference weight sum` = x$reference_weight_sum,
    `Boundary quantile` = x$boundary_quantile,
    `Boundary distance` = x$boundary_distance,
    check.names = FALSE
  )
  print(settings, row.names = FALSE, digits = digits)

  metrics <- rbind(
    `Climatic Displacement` = x$climate_change_amount,
    `Niche Distance Shift` = x$niche_distance_change,
    `Climatic Reconfiguration` = x$climate_reconfiguration,
    `Niche Boundary Exceedance` = x$niche_boundary_exceedance
  )
  cat("\nMetrics\n")
  print(metrics, digits = digits)

  cat("\nNiche Distance Shift direction\n")
  print(x$radial_direction)
  cat("\nNiche boundary status\n")
  print(x$boundary_status)
  invisible(x)
}

#' @export
plot.climniche_fit <- function(x, type = c("distance", "boundary", "amount",
                                           "reconfiguration"), ...) {
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
  }
  invisible(x)
}
