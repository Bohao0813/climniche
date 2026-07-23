#' Extract a tidy climniche table
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells.
#'
#' @return A data frame with one row per evaluated cell, including
#'   `occupied_weight`, the four reported quantities and the old field aliases.
#' @export
climniche_table <- function(x, scope = c("current", "all")) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  idx <- if (scope == "current") x$occupied else seq_len(nrow(x$current))
  occupied_weight <- .fit_reference_weights(x)[idx]
  climate_reconfiguration <- .fit_metric(x, "climate_reconfiguration")[idx]
  niche_boundary_exceedance <- .fit_metric(x, "niche_boundary_exceedance")[idx]
  descriptors <- .fit_exposure_descriptors(x)

  cell <- if (!is.null(x$cell_index) &&
              length(x$cell_index) == nrow(x$current)) {
    x$cell_index[idx]
  } else {
    idx
  }

  data.frame(
    cell = cell,
    occupied_weight = occupied_weight,
    current_niche_distance = x$niche_radius_current[idx],
    future_niche_distance = x$niche_radius_future[idx],
    climate_change_amount = x$climate_change_amount[idx],
    niche_distance_change = x$niche_distance_change[idx],
    climate_reconfiguration = climate_reconfiguration,
    composition_change = climate_reconfiguration,
    change_alignment = x$change_alignment[idx],
    niche_boundary_exceedance = niche_boundary_exceedance,
    outside_niche_exceedance = niche_boundary_exceedance,
    current_niche_percentile = x$niche_percentile$current[idx],
    future_niche_percentile = x$niche_percentile$future[idx],
    percentile_change = x$niche_percentile$delta[idx],
    radial_direction = descriptors$radial_direction[idx],
    boundary_status = descriptors$boundary_status[idx],
    stringsAsFactors = FALSE
  )
}

.summary_weights <- function(tab, scope) {
  if (identical(scope, "current")) {
    tab$occupied_weight
  } else {
    rep(1, nrow(tab))
  }
}

.weighted_prop <- function(x, weights) {
  ok <- !is.na(x) & weights > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(weights[ok] * x[ok]) / sum(weights[ok])
}

.weighted_class_prop <- function(class, weights) {
  lev <- levels(class)
  out <- stats::setNames(rep(0, length(lev)), lev)
  ok <- !is.na(class) & weights > 0
  if (!any(ok)) {
    return(out)
  }
  totals <- tapply(weights[ok], droplevels(class[ok]), sum)
  out[names(totals)] <- totals / sum(weights[ok])
  out[is.na(out)] <- 0
  out
}

.descriptor_summary <- function(tab, weights) {
  radial <- .weighted_class_prop(tab$radial_direction, weights)
  boundary <- .weighted_class_prop(tab$boundary_status, weights)
  data.frame(
    descriptor = c(
      rep("Niche Distance Shift direction", length(radial)),
      rep("Niche boundary status", length(boundary))
    ),
    level = c(names(radial), names(boundary)),
    proportion = c(as.numeric(radial), as.numeric(boundary)),
    stringsAsFactors = FALSE
  )
}

.descriptor_settings_df <- function(x) {
  settings <- x$descriptor_settings %||% x$threshold_settings
  if (is.null(settings)) {
    settings <- list()
  }
  data.frame(
    tolerance = settings$tolerance %||% NA_real_,
    tolerance_quantile = settings$tolerance_quantile %||% NA_real_,
    boundary_exceedance_tolerance =
      settings$boundary_exceedance_tolerance %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.metric_definitions <- function() {
  data.frame(
    metric = c(
      "Climatic Displacement",
      "Niche Distance Shift",
      "Climatic Reconfiguration",
      "Niche Boundary Exceedance"
    ),
    field = c(
      "climate_change_amount",
      "niche_distance_change",
      "climate_reconfiguration",
      "niche_boundary_exceedance"
    ),
    definition = c(
      paste(
        "Distance between current and future conditions under the fitted",
        "climatic metric."
      ),
      paste(
        "Signed change in distance from the current realised climatic",
        "niche centre."
      ),
      paste(
        "Non radial component of Climatic Displacement not captured by",
        "Niche Distance Shift."
      ),
      paste(
        "Positive excess of future niche distance beyond the empirical",
        "radial boundary of the current realised climatic niche."
      )
    ),
    stringsAsFactors = FALSE
  )
}

.report_metric_summary <- function(summary) {
  data.frame(
    metric = c(
      "Climatic Displacement",
      "Climatic Displacement",
      "Climatic Displacement",
      "Niche Distance Shift",
      "Niche Distance Shift",
      "Niche Distance Shift",
      "Climatic Reconfiguration",
      "Climatic Reconfiguration",
      "Climatic Reconfiguration",
      "Niche Boundary Exceedance",
      "Niche Boundary Exceedance",
      "Niche Boundary Exceedance"
    ),
    statistic = c("mean", "median", "90th percentile",
                  "mean", "median", "90th percentile",
                  "mean", "median", "90th percentile",
                  "mean", "median", "90th percentile"),
    value = c(
      summary$mean_climate_change_amount,
      summary$median_climate_change_amount,
      summary$q90_climate_change_amount,
      summary$mean_niche_distance_change,
      summary$median_niche_distance_change,
      summary$q90_niche_distance_change,
      summary$mean_climate_reconfiguration,
      summary$median_climate_reconfiguration,
      summary$q90_climate_reconfiguration,
      summary$mean_niche_boundary_exceedance,
      summary$median_niche_boundary_exceedance,
      summary$q90_niche_boundary_exceedance
    ),
    stringsAsFactors = FALSE
  )
}

.report_variable_display <- function(x) {
  out <- x[, c(
    "variable",
    "mean_absolute_share",
    "mean_signed_contribution",
    "dominant_weight_fraction",
    "interpretation"
  ), drop = FALSE]
  names(out) <- c(
    "Climate variable",
    "Mean absolute share",
    "Mean signed contribution",
    "Weighted dominance frequency",
    "Interpretation"
  )
  out
}

#' Summarise climniche results
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells. Current-scope summaries use reference weights.
#'
#' @return A one-row data frame with summaries of the four continuous reported
#'   quantities, descriptor proportions and fitted descriptor thresholds.
#' @export
climniche_summary <- function(x, scope = c("current", "all")) {
  scope <- match.arg(scope)
  tab <- climniche_table(x, scope = scope)
  weights <- .summary_weights(tab, scope)
  descriptor_summary <- .descriptor_summary(tab, weights)
  get_descriptor_prop <- function(label) {
    value <- descriptor_summary$proportion[descriptor_summary$level == label]
    if (length(value)) value[[1]] else 0
  }
  settings <- .descriptor_settings_df(x)
  boundary_tol <- settings$boundary_exceedance_tolerance
  if (!is.finite(boundary_tol)) {
    boundary_tol <- 0
  }

  out <- cbind(data.frame(
    scope = scope,
    n = nrow(tab),
    reference_weight_sum = sum(weights, na.rm = TRUE),
    boundary_quantile = x$boundary_quantile,
    boundary_distance = x$boundary_radius,
    mean_climate_change_amount =
      .weighted_mean_vector(tab$climate_change_amount, weights),
    median_climate_change_amount = .weighted_quantile(
      tab$climate_change_amount, weights, 0.50
    ),
    q90_climate_change_amount = .weighted_quantile(
      tab$climate_change_amount, weights, 0.90
    ),
    mean_niche_distance_change =
      .weighted_mean_vector(tab$niche_distance_change, weights),
    median_niche_distance_change = .weighted_quantile(
      tab$niche_distance_change, weights, 0.50
    ),
    q90_niche_distance_change = .weighted_quantile(
      tab$niche_distance_change, weights, 0.90
    ),
    mean_climate_reconfiguration =
      .weighted_mean_vector(tab$climate_reconfiguration, weights),
    median_climate_reconfiguration = .weighted_quantile(
      tab$climate_reconfiguration, weights, 0.50
    ),
    q90_climate_reconfiguration = .weighted_quantile(
      tab$climate_reconfiguration, weights, 0.90
    ),
    prop_toward_niche_center = get_descriptor_prop(
      "Toward realised niche centre"
    ),
    prop_limited_niche_distance_shift = get_descriptor_prop(
      "Limited Niche Distance Shift"
    ),
    prop_away_from_niche_center = get_descriptor_prop(
      "Away from realised niche centre"
    ),
    prop_beyond_niche_boundary = get_descriptor_prop(
      "Beyond empirical niche boundary"
    ),
    prop_outside_niche = .weighted_prop(
      tab$niche_boundary_exceedance > boundary_tol, weights
    ),
    mean_niche_boundary_exceedance =
      .weighted_mean_vector(tab$niche_boundary_exceedance, weights),
    median_niche_boundary_exceedance = .weighted_quantile(
      tab$niche_boundary_exceedance, weights, 0.50
    ),
    q90_niche_boundary_exceedance = .weighted_quantile(
      tab$niche_boundary_exceedance, weights, 0.90
    ),
    mean_outside_niche_exceedance =
      .weighted_mean_vector(tab$niche_boundary_exceedance, weights),
    stringsAsFactors = FALSE
  ), settings)
  class(out) <- c("climniche_summary", "data.frame")
  out
}

#' Print a climniche summary
#'
#' Prints the four continuous reported quantities.
#'
#' @param x An object returned by [climniche_summary()].
#' @param ... Additional arguments passed to the data-frame print method.
#'
#' @return Invisibly returns `x`.
#' @export
print.climniche_summary <- function(x, ...) {
  settings <- data.frame(
    Scope = x$scope,
    Cells = x$n,
    `Reference weight sum` = x$reference_weight_sum,
    `Boundary quantile` = x$boundary_quantile,
    `Boundary distance` = x$boundary_distance,
    check.names = FALSE
  )
  metrics <- data.frame(
    Metric = c(
      "Climatic Displacement",
      "Niche Distance Shift",
      "Climatic Reconfiguration",
      "Niche Boundary Exceedance"
    ),
    Mean = c(
      x$mean_climate_change_amount,
      x$mean_niche_distance_change,
      x$mean_climate_reconfiguration,
      x$mean_niche_boundary_exceedance
    ),
    Median = c(
      x$median_climate_change_amount,
      x$median_niche_distance_change,
      x$median_climate_reconfiguration,
      x$median_niche_boundary_exceedance
    ),
    `90th percentile` = c(
      x$q90_climate_change_amount,
      x$q90_niche_distance_change,
      x$q90_climate_reconfiguration,
      x$q90_niche_boundary_exceedance
    ),
    check.names = FALSE
  )
  print.data.frame(settings, row.names = FALSE, ...)
  cat("\n")
  print.data.frame(metrics, row.names = FALSE, ...)
  invisible(x)
}

#' Build a climniche report
#'
#' @param x A fitted `climniche_fit` object.
#' @param species Optional species name used in printed reports.
#' @param scope `"current"` for current reference cells or `"all"` for all
#'   evaluated cells. Current-scope summaries use reference weights.
#' @param top_variables Number of variable contributions to show.
#'
#' @return An object of class `climniche_report`.
#' @export
climniche_report <- function(x, species = NULL, scope = c("current", "all"),
                             top_variables = 5) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  top_variables <- .check_positive_integer(top_variables, "top_variables")
  tab <- climniche_table(x, scope = scope)
  summ <- climniche_summary(x, scope = scope)
  idx <- if (scope == "current") x$occupied else seq_len(nrow(x$current))
  weights <- if (scope == "current") {
    .fit_reference_weights(x)[idx]
  } else {
    rep(1, length(idx))
  }

  descriptor_summary <- .descriptor_summary(tab, weights)

  contribution_fit <- x
  contribution_fit$rasters <- NULL
  contribution_fit$raster_complete <- NULL
  contribution_summary <- climniche_dominant_contribution(
    contribution_fit,
    scope = scope
  )$summary
  var_tab <- data.frame(
    variable = contribution_summary$variable,
    mean_absolute_share = contribution_summary$mean_absolute_share,
    mean_signed_contribution =
      contribution_summary$mean_signed_contribution,
    dominant_weight_fraction =
      contribution_summary$dominant_weight_fraction,
    mean_contribution = contribution_summary$mean_signed_contribution,
    abs_mean_contribution = abs(
      contribution_summary$mean_signed_contribution
    ),
    interpretation = ifelse(
      contribution_summary$mean_signed_contribution > 0,
      "positive mean contribution",
      ifelse(
        contribution_summary$mean_signed_contribution < 0,
        "negative mean contribution",
        "zero mean contribution"
      )
    ),
    stringsAsFactors = FALSE
  )
  var_tab <- var_tab[order(var_tab$mean_absolute_share, decreasing = TRUE),
                     , drop = FALSE]
  var_tab <- utils::head(var_tab, top_variables)

  direction <- if (summ$mean_niche_distance_change > 0) {
    paste(
      "positive, indicating average movement farther from the current",
      "realised climatic niche centre"
    )
  } else if (summ$mean_niche_distance_change < 0) {
    paste(
      "negative, indicating average movement closer to the current",
      "realised climatic niche centre"
    )
  } else {
    "near zero on average"
  }

  exceedance_unit <- if (scope == "current") {
    "total reference weight"
  } else {
    "analysed cells"
  }
  interpretation <- c(
    paste0("Mean Niche Distance Shift is ", direction, "."),
    paste0("Niche Boundary Exceedance is above the configured tolerance for ",
           round(100 * summ$prop_beyond_niche_boundary, 1), "% of ",
           exceedance_unit, " (boundary quantile q = ",
           summ$boundary_quantile, ").")
  )
  metric_weights <- if (!is.null(x$sensitivity_weights)) {
    data.frame(
      variable = names(x$sensitivity_weights),
      sensitivity_weight = as.numeric(x$sensitivity_weights),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      variable = character(),
      sensitivity_weight = numeric(),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    species = species,
    scope = scope,
    metric_definitions = .metric_definitions(),
    settings = data.frame(
      metric_type = x$metric_type %||% "unspecified",
      scaled = x$standardization$enabled %||% NA,
      preprocessed = x$preprocessing$settings$enabled %||% NA,
      boundary_quantile = x$boundary_quantile,
      boundary_distance = x$boundary_distance %||%
        x$boundary_radius %||% x$boundary_value,
      n_variables = ncol(x$current),
      n_cells = nrow(tab),
      stringsAsFactors = FALSE
    ),
    descriptor_settings = .descriptor_settings_df(x),
    summary = summ,
    metric_summary = .report_metric_summary(summ),
    descriptor_summary = descriptor_summary,
    metric_weights = metric_weights,
    metric_matrix = x$A,
    top_variables = var_tab,
    interpretation = interpretation,
    table = tab
  )
  class(out) <- "climniche_report"
  out
}

#' @export
print.climniche_report <- function(x, ...) {
  title <- if (is.null(x$species)) {
    "climniche report"
  } else {
    paste0("climniche report: ", x$species)
  }
  cat(title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")
  cat("Scope: ", x$scope, "\n", sep = "")
  cat("Cells: ", x$settings$n_cells, "\n", sep = "")
  cat("Boundary quantile: ", x$settings$boundary_quantile, "\n\n", sep = "")
  cat("Descriptor settings\n")
  print(x$descriptor_settings, row.names = FALSE)
  cat("\n")

  cat("Interpretation\n")
  for (line in x$interpretation) {
    cat("- ", line, "\n", sep = "")
  }

  cat("\nMetric summary\n")
  print(x$metric_summary, row.names = FALSE)

  if (nrow(x$metric_weights)) {
    cat("\nSensitivity weights\n")
    print(x$metric_weights, row.names = FALSE)
  }

  cat("\nTop variable contributions\n")
  print(.report_variable_display(x$top_variables), row.names = FALSE)
  invisible(x)
}

#' Write a climniche report to Markdown
#'
#' @param report An object returned by [climniche_report()].
#' @param file Output Markdown file.
#'
#' @return Invisibly returns `file`.
#' @export
write_climniche_report <- function(report, file) {
  if (!inherits(report, "climniche_report")) {
    stop("report must be produced by climniche_report().", call. = FALSE)
  }
  title <- if (is.null(report$species)) {
    "Niche relative climate exposure report"
  } else {
    paste0("Niche relative climate exposure report: ",
           report$species)
  }

  fmt_row <- function(dat) {
    paste(utils::capture.output(print(dat, row.names = FALSE)), collapse = "\n")
  }

  lines <- c(
    paste0("# ", title),
    "",
    "## Interpretation",
    paste0("- ", report$interpretation),
    "",
    "## Metric definitions",
    "```text",
    fmt_row(report$metric_definitions),
    "```",
    "",
    "## Fitted settings",
    "```text",
    fmt_row(report$settings),
    "```",
    "",
    "## Descriptor settings",
    "```text",
    fmt_row(report$descriptor_settings),
    "```",
    "",
    "## Metric summary",
    "```text",
    fmt_row(report$metric_summary),
    "```",
    "",
    if (nrow(report$metric_weights)) "## Sensitivity weights" else NULL,
    if (nrow(report$metric_weights)) "```text" else NULL,
    if (nrow(report$metric_weights)) fmt_row(report$metric_weights) else NULL,
    if (nrow(report$metric_weights)) "```" else NULL,
    if (nrow(report$metric_weights)) "" else NULL,
    "## Top Variable Contributions",
    "```text",
    fmt_row(.report_variable_display(report$top_variables)),
    "```",
    "",
    "## Notes",
    paste(
      "The report describes climatic exposure relative to the current",
      "realised climatic niche represented by the supplied occurrence, range",
      "or suitability weights."
    )
  )
  writeLines(lines, con = file)
  invisible(file)
}
