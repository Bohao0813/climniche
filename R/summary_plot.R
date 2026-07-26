#' Summarise a climniche fit
#'
#' @param object A fitted `climniche_fit` object.
#' @param scope `"current"` for a reference-weighted summary or `"all"` for an
#'   unweighted summary across all evaluated rows.
#' @param ... Unused.
#'
#' @return A `summary.climniche_fit` object.
#'
#' @examples
#' sim <- simulate_climniche(n = 200, p = 6, seed = 5)
#' fit <- fit_climniche(
#'   sim$current,
#'   sim$future_away,
#'   occupied = sim$occupied,
#'   sensitivity = sim$sensitivity
#' )
#' summary(fit)
#' @export
summary.climniche_fit <- function(object, scope = c("current", "all"), ...) {
  scope <- match.arg(scope)
  descriptors <- .fit_exposure_descriptors(object)
  reference_weight <- .fit_reference_weights(object)
  summary_weight <- if (identical(scope, "current")) {
    reference_weight
  } else {
    NULL
  }
  summarised <- if (identical(scope, "current")) {
    reference_weight > 0
  } else {
    rep(TRUE, nrow(object$current))
  }
  boundary_distance <- object$boundary_distance %||%
    object$boundary_radius %||% object$boundary_value
  out <- list(
    n = sum(summarised),
    analysed_n = nrow(object$current),
    p = ncol(object$current),
    occupied_n = length(object$occupied),
    reference_weight_sum = sum(reference_weight),
    scope = scope,
    metric_type = object$metric_type %||% "unspecified",
    boundary_quantile = object$boundary_quantile,
    boundary_distance = boundary_distance,
    boundary_value = boundary_distance,
    climate_change_amount = .climniche_numeric_summary(
      object$climate_change_amount, summary_weight
    ),
    niche_distance_change = .climniche_numeric_summary(
      object$niche_distance_change, summary_weight
    ),
    climate_reconfiguration =
      .climniche_numeric_summary(
        .fit_metric(object, "climate_reconfiguration"), summary_weight
      ),
    composition_change = .climniche_numeric_summary(
      .fit_metric(object, "composition_change"), summary_weight
    ),
    change_alignment = .climniche_numeric_summary(
      object$change_alignment, summary_weight
    ),
    niche_boundary_exceedance =
      .climniche_numeric_summary(
        .fit_metric(object, "niche_boundary_exceedance"), summary_weight
      ),
    outside_niche_exceedance =
      .climniche_numeric_summary(
        .fit_metric(object, "outside_niche_exceedance"), summary_weight
      ),
    descriptor_settings = object$descriptor_settings,
    radial_direction = .climniche_descriptor_summary(
      descriptors$radial_direction, summary_weight
    ),
    boundary_status = .climniche_descriptor_summary(
      descriptors$boundary_status, summary_weight
    )
  )
  class(out) <- "summary.climniche_fit"
  out
}

.climniche_numeric_summary <- function(values, weights = NULL) {
  if (is.null(weights)) {
    return(base::summary(values))
  }
  ok <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    out <- stats::setNames(rep(NA_real_, 6L), c(
      "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."
    ))
  } else {
    quantiles <- .weighted_quantile(
      values[ok],
      weights[ok],
      probs = c(0, 0.25, 0.5, 0.75, 1),
      names = FALSE
    )
    out <- c(
      "Min." = quantiles[1L],
      "1st Qu." = quantiles[2L],
      "Median" = quantiles[3L],
      "Mean" = sum(values[ok] * weights[ok]) / sum(weights[ok]),
      "3rd Qu." = quantiles[4L],
      "Max." = quantiles[5L]
    )
  }
  class(out) <- c("summaryDefault", "table")
  out
}

.climniche_descriptor_summary <- function(values, weights = NULL) {
  if (is.null(weights)) {
    return(table(values))
  }
  vapply(levels(values), function(level) {
    sum(weights[!is.na(values) & values == level], na.rm = TRUE)
  }, numeric(1))
}

#' @export
print.summary.climniche_fit <- function(x, digits = max(
  3L, getOption("digits") - 3L
), ...) {
  cat("climniche fit summary\n\n")
  cat(
    "Scope:",
    if (identical(x$scope, "current")) {
      "current reference (weighted)"
    } else {
      "all evaluated rows (unweighted)"
    },
    "\n\n"
  )
  settings <- data.frame(
    `Analysed rows` = x$analysed_n %||% x$n,
    `Summarised rows` = x$n,
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
  cat("\nReported quantities\n")
  print(metrics, digits = digits)

  cat("\nNiche Distance Shift direction\n")
  print(x$radial_direction)
  cat("\nNiche boundary status\n")
  print(x$boundary_status)
  invisible(x)
}

#' Plot a fitted climniche quantity
#'
#' @param x A fitted `climniche_fit` object.
#' @param type Quantity to draw.
#' @param scope `"current"` for a reference-weighted histogram or `"all"` for an
#'   unweighted histogram across evaluated rows.
#' @param breaks Histogram breaks passed to [graphics::hist()].
#' @param ... Additional arguments passed to the histogram plot.
#'
#' @return Invisibly returns `x`.
#' @export
plot.climniche_fit <- function(x, type = c("distance", "boundary", "amount",
                                           "reconfiguration"),
                               scope = c("current", "all"),
                               breaks = "Sturges", ...) {
  type <- match.arg(type)
  scope <- match.arg(scope)
  fields <- c(
    distance = "niche_distance_change",
    boundary = "niche_boundary_exceedance",
    amount = "climate_change_amount",
    reconfiguration = "climate_reconfiguration"
  )
  labels <- c(
    distance = "Niche Distance Shift",
    boundary = "Niche Boundary Exceedance",
    amount = "Climatic Displacement",
    reconfiguration = "Climatic Reconfiguration"
  )
  values <- .fit_metric(x, unname(fields[[type]]))
  weights <- if (identical(scope, "current")) {
    .fit_reference_weights(x)
  } else {
    rep(1, length(values))
  }
  .plot_weighted_histogram(
    values = values,
    weights = weights,
    breaks = breaks,
    main = unname(labels[[type]]),
    xlab = unname(labels[[type]]),
    ylab = if (identical(scope, "current")) {
      "Weighted count"
    } else {
      "Number of cells"
    },
    ...
  )
  invisible(x)
}

.plot_weighted_histogram <- function(values, weights, breaks, main, xlab,
                                     ylab, ...) {
  ok <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    stop("No finite values remain in the selected scope.", call. = FALSE)
  }
  values <- values[ok]
  weights <- weights[ok]
  histogram <- graphics::hist(values, breaks = breaks, plot = FALSE)
  bin <- cut(
    values,
    breaks = histogram$breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  counts <- numeric(length(histogram$counts))
  for (i in seq_along(counts)) {
    counts[i] <- sum(weights[bin == i])
  }
  histogram$counts <- counts
  histogram$density <- counts / sum(counts) / diff(histogram$breaks)
  histogram$intensities <- histogram$density
  graphics::plot(
    histogram,
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  invisible(histogram)
}
