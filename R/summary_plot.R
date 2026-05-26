#' @export
summary.climniche_fit <- function(object, ...) {
  out <- list(
    n = nrow(object$current),
    p = ncol(object$current),
    occupied_n = length(object$occupied),
    boundary_quantile = object$boundary_quantile,
    boundary_value = object$boundary_value,
    climate_change_amount = base::summary(object$climate_change_amount),
    niche_distance_change = base::summary(object$niche_distance_change),
    composition_change = base::summary(object$composition_change),
    change_alignment = base::summary(object$change_alignment),
    outside_niche_exceedance = base::summary(object$outside_niche_exceedance),
    classification = table(object$classification)
  )
  class(out) <- "summary.climniche_fit"
  out
}

#' @export
plot.climniche_fit <- function(x, type = c("distance", "boundary", "amount",
                                           "classification"), ...) {
  type <- match.arg(type)
  if (type == "distance") {
    graphics::hist(x$niche_distance_change, main = "Niche distance change",
                   xlab = "Niche distance change", ...)
  } else if (type == "boundary") {
    graphics::hist(x$outside_niche_exceedance, main = "Niche boundary exceedance",
                   xlab = "Exceedance", ...)
  } else if (type == "amount") {
    graphics::hist(x$climate_change_amount, main = "Climate change amount",
                   xlab = "Amount", ...)
  } else {
    graphics::barplot(table(x$classification), las = 2,
                      main = "Climate niche change classes", ...)
  }
  invisible(x)
}
