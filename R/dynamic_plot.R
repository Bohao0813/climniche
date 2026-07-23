.series_plot_labels <- function(metric) {
  switch(
    metric,
    exposed_fraction = "Weighted Fraction Beyond Niche Boundary",
    conditional_relative_exceedance = "Conditional Boundary Exceedance",
    range_wide_relative_exceedance =
      "Range Mean Relative Boundary Exceedance",
    mean_niche_boundary_exceedance = "Mean Niche Boundary Exceedance",
    metric
  )
}

.series_palette <- function(n) {
  base <- c("#246a73", "#b64f4a", "#d28e2f", "#3f7d55",
            "#725a9b", "#3e77a8")
  rep(base, length.out = n)
}

.series_plot_summary <- function(data, metric, interval) {
  scenario <- ifelse(is.na(data$scenario), "Projection", data$scenario)
  groups <- split(
    seq_len(nrow(data)),
    paste(scenario, as.character(data$time), sep = "\r")
  )
  alpha <- (1 - interval) / 2
  rows <- lapply(groups, function(index) {
    dat <- data[index, , drop = FALSE]
    values <- dat[[metric]]
    values <- values[is.finite(values)]
    q <- if (length(values)) {
      stats::quantile(
        values,
        probs = c(alpha, 0.5, 1 - alpha),
        names = FALSE,
        type = 8
      )
    } else {
      rep(NA_real_, 3L)
    }
    data.frame(
      time = dat$time[1L],
      scenario_plot = ifelse(
        is.na(dat$scenario[1L]),
        "Projection",
        dat$scenario[1L]
      ),
      lower = q[1L],
      median = q[2L],
      upper = q[3L],
      n_models = length(values),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Plot range wide niche boundary exceedance through time
#'
#' @param x A `climniche_series` object.
#' @param metric Range summary to display.
#' @param scope,aggregation_weight,area_weight,boundary_exceedance_tolerance
#'   Arguments passed to [climniche_range_summary()].
#' @param interval Central interval drawn across climate models.
#' @param show_models If `TRUE`, draw individual climate model lines behind the
#'   ensemble median.
#' @param title Optional plot title.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @export
plot_climniche_time <- function(
    x,
    metric = c("range_wide_relative_exceedance",
               "exposed_fraction",
               "conditional_relative_exceedance",
               "mean_niche_boundary_exceedance"),
    scope = c("current", "all"),
    aggregation_weight = NULL,
    area_weight = FALSE,
    boundary_exceedance_tolerance = NULL,
    interval = 0.80,
    show_models = TRUE,
    title = NULL,
    base_size = 9) {
  .need_ggplot2()
  if (!inherits(x, "climniche_series")) {
    stop("x must be a climniche_series object.", call. = FALSE)
  }
  metric <- match.arg(metric)
  scope <- match.arg(scope)
  interval <- .check_open_probability(interval, "interval")
  show_models <- .check_flag(show_models, "show_models")
  base_size <- .check_finite_scalar(base_size, "base_size")
  if (base_size <= 0) {
    stop("base_size must be positive.", call. = FALSE)
  }
  data <- climniche_range_summary(
    x,
    scope = scope,
    aggregation_weight = aggregation_weight,
    area_weight = area_weight,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance
  )
  data$scenario_plot <- ifelse(
    is.na(data$scenario),
    "Projection",
    data$scenario
  )
  data$model_plot <- ifelse(is.na(data$model), "Model", data$model)
  data$plot_value <- data[[metric]]
  ensemble <- .series_plot_summary(data, metric, interval)
  scenarios <- unique(ensemble$scenario_plot)
  colours <- stats::setNames(.series_palette(length(scenarios)), scenarios)
  multiple_scenarios <- length(scenarios) > 1L ||
    !identical(scenarios, "Projection")
  multiple_models <- any(ensemble$n_models > 1L)

  p <- ggplot2::ggplot()
  if (show_models && multiple_models) {
    p <- p + ggplot2::geom_line(
      data = data,
      ggplot2::aes(
        x = time,
        y = plot_value,
        colour = scenario_plot,
        group = interaction(scenario_plot, model_plot, drop = TRUE)
      ),
      linewidth = 0.35,
      alpha = 0.28
    )
  }
  if (multiple_models) {
    p <- p + ggplot2::geom_ribbon(
      data = ensemble,
      ggplot2::aes(
        x = time,
        ymin = lower,
        ymax = upper,
        fill = scenario_plot,
        group = scenario_plot
      ),
      alpha = 0.18,
      colour = NA
    )
  }
  p <- p + ggplot2::geom_line(
    data = ensemble,
    ggplot2::aes(
      x = time,
      y = median,
      colour = scenario_plot,
      group = scenario_plot
    ),
    linewidth = 0.8
  ) +
    ggplot2::geom_point(
      data = ensemble,
      ggplot2::aes(x = time, y = median, colour = scenario_plot),
      size = 1.35
    ) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(
      x = "Projection Time",
      y = .series_plot_labels(metric),
      colour = if (multiple_scenarios) "Scenario" else NULL,
      title = title
    ) +
    .climniche_theme(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (multiple_scenarios) "right" else "none",
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey88",
        linewidth = 0.25
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
  if (multiple_models) {
    p <- p +
      ggplot2::scale_fill_manual(values = colours) +
      ggplot2::labs(
        fill = if (multiple_scenarios) "Scenario" else NULL
      )
  }
  if (identical(metric, "exposed_fraction")) {
    p <- p + ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = function(z) paste0(round(100 * z), "%")
    )
  }
  p
}

.select_series_value <- function(values, requested, arg) {
  available <- unique(values[!is.na(values)])
  if (!is.null(requested)) {
    requested <- as.character(requested)[1L]
    if (!requested %in% available) {
      stop(arg, " was not found in the fitted series.", call. = FALSE)
    }
    return(requested)
  }
  if (length(available) > 1L) {
    stop("Select one ", arg, " for this map.", call. = FALSE)
  }
  if (length(available)) available else NA_character_
}

.select_series_time <- function(values, requested) {
  available <- unique(values)
  if (!is.null(requested)) {
    match_index <- match(as.character(requested), as.character(available))
    if (is.na(match_index)) {
      stop("time was not found in the fitted series.", call. = FALSE)
    }
    return(available[match_index])
  }
  if (length(available) > 1L) {
    stop("Select one time for this map.", call. = FALSE)
  }
  available[1L]
}

.series_filter <- function(data, model = NULL, scenario = NULL,
                           time = NULL, use_model = TRUE) {
  selected_scenario <- .select_series_value(
    data$scenario,
    scenario,
    "scenario"
  )
  keep <- if (is.na(selected_scenario)) {
    is.na(data$scenario)
  } else {
    data$scenario == selected_scenario
  }
  data <- data[keep, , drop = FALSE]
  if (use_model) {
    selected_model <- .select_series_value(data$model, model, "model")
    keep <- if (is.na(selected_model)) {
      is.na(data$model)
    } else {
      data$model == selected_model
    }
    data <- data[keep, , drop = FALSE]
  }
  if (!is.null(time) || length(unique(data$time)) > 1L) {
    selected_time <- .select_series_time(data$time, time)
    data <- data[as.character(data$time) == as.character(selected_time),
                 , drop = FALSE]
  }
  data
}

#' Map dynamic niche boundary departure
#'
#' @param x A spatial `climniche_series` object.
#' @param metric Dynamic result to map.
#' @param model,scenario,time Optional projection selectors. A selector is
#'   required when more than one relevant value is present.
#' @param scope,persistence,boundary_exceedance_tolerance Arguments passed to
#'   [climniche_departure()] or [climniche_model_agreement()].
#' @param title Optional map title. Use `FALSE` to suppress it.
#' @param legend_title Optional legend title. The default suppresses a title
#'   that would repeat the map title.
#' @param limits,colours Optional colour-scale limits and colours.
#' @param ... Additional arguments passed to [plot_climniche_map()].
#'
#' @return A ggplot object.
#' @export
plot_climniche_departure_map <- function(
    x,
    metric = c("first_persistent_departure",
               "departure_time_fraction",
               "mean_relative_exceedance",
               "maximum_relative_exceedance",
               "model_agreement"),
    model = NULL,
    scenario = NULL,
    time = NULL,
    scope = c("current", "all"),
    persistence = 1L,
    boundary_exceedance_tolerance = NULL,
    title = NULL,
    legend_title = FALSE,
    limits = NULL,
    colours = NULL,
    ...) {
  .need_ggplot2()
  if (!inherits(x, "climniche_series") ||
      identical(x$input_type, "matrix")) {
    stop("x must be a spatial climniche_series object.", call. = FALSE)
  }
  metric <- match.arg(metric)
  scope <- match.arg(scope)
  fit <- x$fits[[1L]]

  if (identical(metric, "model_agreement")) {
    data <- climniche_model_agreement(
      x,
      scope = scope,
      boundary_exceedance_tolerance =
        boundary_exceedance_tolerance
    )
    data <- .series_filter(
      data,
      scenario = scenario,
      time = time,
      use_model = FALSE
    )
    value <- data$model_agreement
    default_title <- "Climate Model Agreement"
    default_legend <- "Model Agreement"
    if (is.null(limits)) limits <- c(0, 1)
    if (is.null(colours)) {
      colours <- c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5")
    }
  } else {
    if (!is.null(time)) {
      stop("time is only used for model_agreement maps.", call. = FALSE)
    }
    data <- climniche_departure(
      x,
      scope = scope,
      persistence = persistence,
      boundary_exceedance_tolerance =
        boundary_exceedance_tolerance
    )
    data <- .series_filter(
      data,
      model = model,
      scenario = scenario,
      use_model = TRUE
    )
    value <- data[[metric]]
    labels <- c(
      first_persistent_departure = "First Projected Persistent Departure",
      departure_time_fraction = "Time Beyond Niche Boundary",
      mean_relative_exceedance = "Mean Relative Boundary Exceedance",
      maximum_relative_exceedance = "Maximum Relative Boundary Exceedance"
    )
    default_title <- labels[[metric]]
    default_legend <- labels[[metric]]
    if (identical(metric, "first_persistent_departure")) {
      if (!is.numeric(x$index$time)) {
        stop("First-departure maps currently require numeric projection times.",
             call. = FALSE)
      }
      value <- as.numeric(value)
    }
    if (identical(metric, "departure_time_fraction") && is.null(limits)) {
      limits <- c(0, 1)
    }
    if (is.null(colours)) {
      colours <- c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5")
    }
  }
  if (!nrow(data) || !any(is.finite(value))) {
    stop("No finite dynamic values are available for the selected map.",
         call. = FALSE)
  }

  fit_cells <- if (is.null(fit$cell_index)) {
    seq_len(nrow(fit$current))
  } else {
    fit$cell_index
  }
  projected_values <- rep(NA_real_, nrow(fit$current))
  cell_match <- match(data$cell, fit_cells)
  projected_values[cell_match[!is.na(cell_match)]] <-
    value[!is.na(cell_match)]
  if (.is_raster(.fit_spatial_template(fit))) {
    raster <- .values_to_raster(
      .fit_spatial_template(fit),
      projected_values,
      fit$raster_complete
    )
  } else {
    raster <- .values_to_spatraster(
      .fit_spatial_template(fit),
      projected_values,
      fit$raster_complete
    )
  }
  if (is.null(title)) title <- default_title
  if (is.null(legend_title)) legend_title <- default_legend
  .plot_climniche_map(
    raster,
    metric = "niche_boundary_exceedance",
    title = title,
    legend_title = legend_title,
    limits = limits,
    colours = colours,
    ...
  )
}

#' @export
plot.climniche_series <- function(x, ...) {
  plot_climniche_time(x, ...)
}
