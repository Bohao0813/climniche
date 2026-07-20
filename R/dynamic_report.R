.series_cell_weights <- function(x, scope, aggregation_weight, area_weight) {
  fit <- x$fits[[1L]]
  reference_weight <- .fit_reference_weights(fit)
  modifier <- .aggregation_modifier(fit, aggregation_weight)
  if (area_weight) {
    modifier <- modifier * .spatial_area_weights(fit)
  }
  weight <- if (identical(scope, "current")) {
    reference_weight * modifier
  } else {
    modifier
  }
  cell <- if (is.null(fit$cell_index)) {
    seq_len(nrow(fit$current))
  } else {
    fit$cell_index
  }
  data.frame(cell = cell, aggregation_weight = weight)
}

.dynamic_definitions <- function() {
  data.frame(
    name = c(
      "Weighted Exposed Fraction",
      "Conditional Relative Boundary Exceedance",
      "Range Mean Relative Boundary Exceedance",
      "First Projected Persistent Departure",
      "Boundary Departure Time Fraction",
      "Cumulative Relative Boundary Exceedance",
      "Climate Model Agreement",
      "Maximum Interval Increase Rate"
    ),
    field = c(
      "exposed_fraction",
      "conditional_relative_exceedance",
      "range_wide_relative_exceedance",
      "first_persistent_departure",
      "departure_time_fraction",
      "cumulative_relative_exceedance",
      "model_agreement",
      "maximum_increase_rate"
    ),
    definition = c(
      "Fraction of the summary weight beyond the empirical niche boundary.",
      paste(
        "Mean boundary exceedance among exposed cells, divided by the",
        "fitted boundary distance."
      ),
      paste(
        "Mean relative boundary exceedance across all cells included in",
        "the range summary."
      ),
      "First sampled projection in a run meeting the persistence setting.",
      paste(
        "Fraction of the sampled time interval beyond the empirical niche",
        "boundary."
      ),
      paste(
        "Trapezoidal time integral of relative boundary exceedance in the",
        "reported time units."
      ),
      paste(
        "Proportion of at least two supplied climate models projecting",
        "positive boundary exceedance."
      ),
      paste(
        "Largest positive finite difference in a range summary divided by",
        "its projection interval."
      )
    ),
    stringsAsFactors = FALSE
  )
}

.series_display_table <- function(data) {
  if (is.null(data)) {
    return(NULL)
  }
  if ("metric" %in% names(data)) {
    metric_labels <- c(
      exposed_fraction = "Weighted Exposed Fraction",
      conditional_relative_exceedance =
        "Conditional Relative Boundary Exceedance",
      range_wide_relative_exceedance =
        "Range Mean Relative Boundary Exceedance",
      mean_niche_boundary_exceedance = "Mean Niche Boundary Exceedance"
    )
    matched_metric <- match(data$metric, names(metric_labels))
    replace_metric <- !is.na(matched_metric)
    data$metric[replace_metric] <-
      unname(metric_labels[matched_metric[replace_metric]])
  }
  labels <- c(
    projection = "Projection",
    time = "Projection Time",
    model = "Climate Model",
    scenario = "Scenario",
    scope = "Scope",
    n = "Cells",
    n_cells = "Cells",
    n_models = "Climate Models",
    n_projections = "Projections",
    aggregation_weight_sum = "Weight",
    boundary_distance = "Niche Boundary Distance",
    boundary_quantile = "Niche Boundary Quantile",
    boundary_exceedance_tolerance = "Boundary Tolerance",
    exposed_fraction = "Weighted Exposed Fraction",
    conditional_relative_exceedance =
      "Conditional Relative Boundary Exceedance",
    range_wide_relative_exceedance =
      "Range Mean Relative Boundary Exceedance",
    mean_niche_boundary_exceedance = "Mean Niche Boundary Exceedance",
    proportion_with_persistent_departure =
      "Fraction With Persistent Departure",
    median_first_persistent_departure =
      "Median First Projected Persistent Departure",
    mean_departure_time_fraction = "Mean Boundary Departure Time Fraction",
    mean_relative_exceedance = "Mean Relative Boundary Exceedance",
    proportion_with_reentry = "Fraction With Re-entry",
    metric = "Metric",
    first_time = "First Projection Time",
    last_time = "Last Projection Time",
    total_change = "Total Change",
    maximum_interval_increase = "Maximum Interval Increase",
    maximum_increase_rate = "Maximum Interval Increase Rate",
    interval_start = "Interval Start",
    interval_end = "Interval End",
    time_unit = "Time Unit",
    mean_model_agreement = "Mean Climate Model Agreement",
    persistence = "Persistence",
    tolerance = "Niche Distance Shift Tolerance",
    tolerance_quantile = "Tolerance Quantile",
    area_weight = "Area Weight",
    agreement_interval = "Model Interval",
    n_reference_rows = "Reference Cells",
    positive_reference_rows = "Positive Reference Cells",
    retained_variables = "Retained Climate Variables"
  )
  matched <- match(names(data), names(labels))
  replace <- !is.na(matched)
  names(data)[replace] <- unname(labels[matched[replace]])
  data
}

.weighted_departure_summary <- function(departure, weights) {
  data <- merge(departure, weights, by = "cell", all.x = TRUE, sort = FALSE)
  groups <- split(
    seq_len(nrow(data)),
    .series_group_key(data$model, data$scenario)
  )
  rows <- lapply(groups, function(index) {
    dat <- data[index, , drop = FALSE]
    weight <- .clean_aggregation_weights(
      dat$aggregation_weight,
      "aggregation weights"
    )
    ok <- weight > 0
    weight <- weight[ok]
    dat <- dat[ok, , drop = FALSE]
    first_numeric <- as.numeric(dat$first_persistent_departure)
    first_ok <- is.finite(first_numeric)
    data.frame(
      model = dat$model[1L],
      scenario = dat$scenario[1L],
      n = nrow(dat),
      aggregation_weight_sum = sum(weight),
      proportion_with_persistent_departure =
        .weighted_mean_vector(
          as.numeric(!is.na(dat$first_persistent_departure)),
          weight
        ),
      median_first_persistent_departure = if (any(first_ok)) {
        .restore_time_class(.weighted_quantile(
          first_numeric[first_ok],
          weight[first_ok],
          probs = 0.5
        ), dat$first_persistent_departure)
      } else {
        NA_real_
      },
      mean_departure_time_fraction = .weighted_mean_vector(
        dat$departure_time_fraction,
        weight
      ),
      mean_relative_exceedance = .weighted_mean_vector(
        dat$mean_relative_exceedance,
        weight
      ),
      proportion_with_reentry = .weighted_mean_vector(
        as.numeric(dat$reentered),
        weight
      ),
      time_unit = dat$time_unit[1L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.weighted_agreement_summary <- function(agreement, weights) {
  data <- merge(agreement, weights, by = "cell", all.x = TRUE, sort = FALSE)
  groups <- split(
    seq_len(nrow(data)),
    .series_group_key(
      model = rep(NA_character_, nrow(data)),
      scenario = data$scenario,
      time = data$time
    )
  )
  rows <- lapply(groups, function(index) {
    dat <- data[index, , drop = FALSE]
    data.frame(
      time = dat$time[1L],
      scenario = dat$scenario[1L],
      n_cells = nrow(dat),
      mean_model_agreement = .weighted_mean_vector(
        dat$model_agreement,
        dat$aggregation_weight
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Build a report for a climniche time series
#'
#' @param x A `climniche_series` object.
#' @param species Optional species or study label.
#' @param scope,aggregation_weight,area_weight,boundary_exceedance_tolerance
#'   Range-summary settings.
#' @param persistence Consecutive projection count used to identify persistent
#'   boundary departure.
#' @param agreement_interval Central interval used for model agreement.
#'
#' @return A `climniche_series_report` object.
#' @export
climniche_series_report <- function(
    x,
    species = NULL,
    scope = c("current", "all"),
    aggregation_weight = NULL,
    area_weight = FALSE,
    boundary_exceedance_tolerance = NULL,
    persistence = 1L,
    agreement_interval = 0.80) {
  if (!inherits(x, "climniche_series")) {
    stop("x must be a climniche_series object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  area_weight <- .check_flag(area_weight, "area_weight")
  persistence <- .check_positive_integer(persistence, "persistence")
  agreement_interval <- .check_open_probability(
    agreement_interval,
    "agreement_interval"
  )
  weights <- .series_cell_weights(
    x,
    scope = scope,
    aggregation_weight = aggregation_weight,
    area_weight = area_weight
  )
  departure <- climniche_departure(
    x,
    scope = scope,
    persistence = persistence,
    boundary_exceedance_tolerance =
      boundary_exceedance_tolerance
  )
  agreement <- climniche_model_agreement(
    x,
    scope = scope,
    boundary_exceedance_tolerance =
      boundary_exceedance_tolerance,
    interval = agreement_interval
  )
  change_rate_metrics <- c(
    "exposed_fraction",
    "range_wide_relative_exceedance",
    "mean_niche_boundary_exceedance"
  )
  change_rate <- do.call(rbind, lapply(change_rate_metrics, function(metric) {
    climniche_change_rate(
      x,
      metric = metric,
      scope = scope,
      aggregation_weight = aggregation_weight,
      area_weight = area_weight,
      boundary_exceedance_tolerance =
        boundary_exceedance_tolerance
    )
  }))
  rownames(change_rate) <- NULL
  agreement_summary <- if (any(agreement$n_models >= 2L)) {
    .weighted_agreement_summary(agreement, weights)
  } else {
    NULL
  }
  out <- list(
    species = species,
    call = match.call(),
    reference = data.frame(
      n_reference_rows = length(x$reference$reference_weight),
      positive_reference_rows = sum(x$reference$reference_weight > 0),
      retained_variables = ncol(x$reference$current),
      boundary_quantile = x$reference$boundary_quantile,
      boundary_distance = x$reference$boundary_distance
    ),
    settings = data.frame(
      scope = scope,
      area_weight = area_weight,
      persistence = persistence,
      tolerance = x$descriptor_settings$tolerance,
      tolerance_quantile = x$descriptor_settings$tolerance_quantile,
      boundary_exceedance_tolerance =
        .range_boundary_tolerance(
          x$fits[[1L]],
          boundary_exceedance_tolerance
        ),
      agreement_interval = agreement_interval
    ),
    projections = x$index,
    range_summary = climniche_range_summary(
      x,
      scope = scope,
      aggregation_weight = aggregation_weight,
      area_weight = area_weight,
      boundary_exceedance_tolerance =
        boundary_exceedance_tolerance
    ),
    departure_summary = .weighted_departure_summary(departure, weights),
    change_rate = change_rate,
    model_agreement_summary = agreement_summary,
    metric_definitions = .metric_definitions(),
    dynamic_definitions = .dynamic_definitions()
  )
  class(out) <- "climniche_series_report"
  out
}

#' @export
print.climniche_series_report <- function(x, ...) {
  if (!is.null(x$species) && nzchar(x$species)) {
    cat("Niche climate exposure series:", x$species, "\n\n")
  } else {
    cat("Niche climate exposure series\n\n")
  }
  cat("Reference\n")
  print(.series_display_table(x$reference), row.names = FALSE)
  cat("\nSettings\n")
  print(.series_display_table(x$settings), row.names = FALSE)
  cat("\nRange through time\n")
  print(.series_display_table(x$range_summary), row.names = FALSE)
  cat("\nDynamic boundary departure\n")
  print(.series_display_table(x$departure_summary), row.names = FALSE)
  cat("\nMaximum interval increases\n")
  print(.series_display_table(x$change_rate), row.names = FALSE)
  if (!is.null(x$model_agreement_summary)) {
    cat("\nClimate model agreement\n")
    print(.series_display_table(x$model_agreement_summary), row.names = FALSE)
  }
  invisible(x)
}

#' Write a climniche time series report
#'
#' @param report A `climniche_series_report` object.
#' @param file Output Markdown file.
#'
#' @return The file path, invisibly.
#' @export
write_climniche_series_report <- function(report, file) {
  if (!inherits(report, "climniche_series_report")) {
    stop("report must be a climniche_series_report object.", call. = FALSE)
  }
  file <- as.character(file)[1L]
  if (is.na(file) || !nzchar(file)) {
    stop("file must be a non-empty path.", call. = FALSE)
  }
  format_table <- function(data) {
    output <- utils::capture.output(print(
      .series_display_table(data),
      row.names = FALSE
    ))
    paste(output, collapse = "\n")
  }
  title <- if (!is.null(report$species) && nzchar(report$species)) {
    paste0("Niche climate exposure series: ", report$species)
  } else {
    "Niche climate exposure series"
  }
  lines <- c(
    paste0("# ", title),
    "",
    "## Reference",
    "```text",
    format_table(report$reference),
    "```",
    "",
    "## Settings",
    "```text",
    format_table(report$settings),
    "```",
    "",
    "## Projections",
    "```text",
    format_table(report$projections),
    "```",
    "",
    "## Range Through Time",
    "```text",
    format_table(report$range_summary),
    "```",
    "",
    "## Dynamic Boundary Departure",
    "```text",
    format_table(report$departure_summary),
    "```",
    "",
    "## Interval Change",
    "```text",
    format_table(report$change_rate),
    "```",
    "",
    if (!is.null(report$model_agreement_summary)) {
      "## Climate Model Agreement"
    } else {
      NULL
    },
    if (!is.null(report$model_agreement_summary)) "```text" else NULL,
    if (!is.null(report$model_agreement_summary)) {
      format_table(report$model_agreement_summary)
    } else {
      NULL
    },
    if (!is.null(report$model_agreement_summary)) "```" else NULL,
    if (!is.null(report$model_agreement_summary)) "" else NULL,
    "## Metric Definitions",
    "```text",
    format_table(report$metric_definitions),
    "```",
    "",
    "## Dynamic Summary Definitions",
    "```text",
    format_table(report$dynamic_definitions),
    "```"
  )
  writeLines(lines, con = file, useBytes = TRUE)
  invisible(file)
}
