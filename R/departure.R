#' Summarise niche boundary departure through time
#'
#' @param x A `climniche_series` object.
#' @param scope `"current"` for positive reference weights or `"all"` for every
#'   evaluated cell.
#' @param persistence Minimum number of consecutive projections required to
#'   identify a persistent departure. This affects persistent timing summaries,
#'   not the continuous exceedance values.
#' @param boundary_exceedance_tolerance Optional non-negative boundary
#'   tolerance. The fitted value is used by default.
#'
#' @return A data frame with one row per cell, model and scenario.
#'
#' @details
#' At projection time \eqn{t_k}, let \eqn{e_{ik}} be Niche Boundary Exceedance
#' for cell \eqn{i}, divided by the fitted niche boundary distance after the
#' boundary tolerance is applied. Cumulative relative exceedance is
#' \deqn{J_i = \sum_{k=1}^{K-1} (t_{k+1}-t_k)
#'       \frac{e_{ik}+e_{i,k+1}}{2}.}
#' Mean relative exceedance is \eqn{J_i/(t_K-t_1)}. Departure time fraction
#' applies the same trapezoidal calculation to the binary boundary state.
#' Numeric, Date and POSIXct times may be irregularly spaced. Cumulative values
#' use the reported time unit; mean values are dimensionless.
#'
#' A persistent departure is a run of at least `persistence` sampled
#' projections beyond the boundary. The reported first departure is the first
#' projection in that run, not an estimated event date between projections.
#'
#' @examples
#' sim <- simulate_climniche(n = 160, p = 6, seed = 2)
#' future <- lapply(c(0.2, 0.5, 0.8, 1), function(fraction) {
#'   sim$current + fraction * (sim$future_away - sim$current)
#' })
#' series <- fit_climniche_series(
#'   sim$current,
#'   future,
#'   time = c(2030, 2050, 2070, 2090),
#'   occupied = sim$occupied,
#'   sensitivity = sim$sensitivity
#' )
#' head(climniche_departure(series))
#' @export
climniche_departure <- function(
    x,
    scope = c("current", "all"),
    persistence = 1L,
    boundary_exceedance_tolerance = NULL) {
  if (!inherits(x, "climniche_series")) {
    stop("x must be a climniche_series object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  persistence <- .check_positive_integer(persistence, "persistence")
  tolerance <- .range_boundary_tolerance(
    x$fits[[1L]],
    boundary_exceedance_tolerance
  )
  tab <- climniche_series_table(x, scope = scope)
  boundary <- as.numeric(x$reference$boundary_radius)[1L]
  if (!is.finite(boundary) || boundary < 0) {
    stop("x contains an invalid niche boundary distance.", call. = FALSE)
  }
  groups <- split(
    seq_len(nrow(tab)),
    .series_group_key(tab$model, tab$scenario, tab$cell)
  )
  rows <- lapply(groups, function(index) {
    dat <- tab[index, , drop = FALSE]
    ord <- order(as.numeric(dat$time))
    dat <- dat[ord, , drop = FALSE]
    time_numeric <- .time_numeric(dat$time)
    exceedance <- dat$niche_boundary_exceedance
    beyond <- is.finite(exceedance) & exceedance > tolerance
    persistent <- .persistent_departure_flags(beyond, persistence)
    relative <- if (boundary > 0) {
      ifelse(beyond, exceedance / boundary, 0)
    } else {
      rep(NA_real_, length(exceedance))
    }
    first_index <- which(persistent)[1L]
    last_index <- utils::tail(which(persistent), 1L)
    first_time <- if (length(first_index) && !is.na(first_index)) {
      dat$time[first_index]
    } else {
      dat$time[NA_integer_]
    }
    last_time <- if (length(last_index)) {
      dat$time[last_index]
    } else {
      dat$time[NA_integer_]
    }
    reentry_count <- if (length(persistent) > 1L) {
      sum(diff(as.integer(persistent)) == -1L)
    } else {
      0L
    }
    data.frame(
      cell = dat$cell[1L],
      model = dat$model[1L],
      scenario = dat$scenario[1L],
      scope = scope,
      n_projections = nrow(dat),
      persistence = persistence,
      first_persistent_departure = first_time,
      last_persistent_departure = last_time,
      departure_projection_fraction = mean(beyond),
      departure_time_fraction = .time_mean(as.numeric(beyond), time_numeric),
      persistent_departure_time_fraction =
        .time_mean(as.numeric(persistent), time_numeric),
      mean_relative_exceedance = .time_mean(relative, time_numeric),
      cumulative_relative_exceedance =
        .time_integral(relative, time_numeric),
      maximum_relative_exceedance = if (all(is.na(relative))) {
        NA_real_
      } else {
        max(relative, na.rm = TRUE)
      },
      reentry_count = reentry_count,
      reentered = reentry_count > 0L,
      boundary_distance = boundary,
      boundary_exceedance_tolerance = tolerance,
      time_unit = .time_unit(dat$time),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("climniche_departure", "data.frame")
  out
}

.series_group_key <- function(model, scenario, cell = NULL, time = NULL) {
  parts <- list(
    .series_label_key(model),
    .series_label_key(scenario)
  )
  if (!is.null(time)) {
    parts <- c(parts, list(as.character(time)))
  }
  if (!is.null(cell)) {
    parts <- c(parts, list(as.character(cell)))
  }
  do.call(paste, c(parts, sep = "\r"))
}

.persistent_departure_flags <- function(x, persistence) {
  out <- rep(FALSE, length(x))
  runs <- rle(x)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  keep <- which(runs$values & runs$lengths >= persistence)
  for (i in keep) {
    out[starts[i]:ends[i]] <- TRUE
  }
  out
}

.time_numeric <- function(x) {
  value <- as.numeric(x)
  if (any(!is.finite(value))) {
    stop("Projection times must be finite.", call. = FALSE)
  }
  value
}

.time_unit <- function(x) {
  if (inherits(x, "Date")) {
    return("days")
  }
  if (inherits(x, "POSIXt")) {
    return("seconds")
  }
  "supplied time units"
}

.restore_time_class <- function(value, template) {
  if (inherits(template, "Date")) {
    return(as.Date(value, origin = "1970-01-01"))
  }
  if (inherits(template, "POSIXt")) {
    timezone <- attr(template, "tzone")
    if (is.null(timezone) || !length(timezone) || !nzchar(timezone[1L])) {
      timezone <- "UTC"
    }
    return(as.POSIXct(value, origin = "1970-01-01", tz = timezone[1L]))
  }
  value
}

.time_integral <- function(y, time) {
  ok <- is.finite(y) & is.finite(time)
  y <- y[ok]
  time <- time[ok]
  if (!length(y)) {
    return(NA_real_)
  }
  if (length(y) == 1L) {
    return(0)
  }
  ord <- order(time)
  y <- y[ord]
  time <- time[ord]
  if (anyDuplicated(time)) {
    stop("Projection times must be unique within each model and scenario.",
         call. = FALSE)
  }
  sum(diff(time) * (utils::head(y, -1L) + utils::tail(y, -1L)) / 2)
}

.time_mean <- function(y, time) {
  ok <- is.finite(y) & is.finite(time)
  y <- y[ok]
  time <- time[ok]
  if (!length(y)) {
    return(NA_real_)
  }
  if (length(y) == 1L || diff(range(time)) == 0) {
    return(mean(y))
  }
  .time_integral(y, time) / diff(range(time))
}

#' Summarise agreement among projected climate models
#'
#' @param x A `climniche_series` object.
#' @param scope `"current"` or `"all"`.
#' @param boundary_exceedance_tolerance Optional non-negative boundary
#'   tolerance.
#' @param interval Central interval reported for Niche Boundary Exceedance
#'   across models.
#'
#' @return A data frame with one row per cell, time and scenario.
#'
#' @details
#' `model_agreement` is the proportion of supplied climate models with positive
#' Niche Boundary Exceedance. It is an ensemble agreement measure, not an
#' occurrence probability. Agreement is `NA` when fewer than two finite model
#' projections are available for a cell, time and scenario.
#' For \eqn{M} available models,
#' \deqn{G_{it} = \frac{1}{M}\sum_{m=1}^{M} I(E_{imt} > \tau),}
#' where \eqn{\tau} is `boundary_exceedance_tolerance`.
#' @export
climniche_model_agreement <- function(
    x,
    scope = c("current", "all"),
    boundary_exceedance_tolerance = NULL,
    interval = 0.80) {
  if (!inherits(x, "climniche_series")) {
    stop("x must be a climniche_series object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  interval <- .check_open_probability(interval, "interval")
  tolerance <- .range_boundary_tolerance(
    x$fits[[1L]],
    boundary_exceedance_tolerance
  )
  tab <- climniche_series_table(x, scope = scope)
  groups <- split(
    seq_len(nrow(tab)),
    .series_group_key(
      model = rep(NA_character_, nrow(tab)),
      scenario = tab$scenario,
      cell = tab$cell,
      time = tab$time
    )
  )
  alpha <- (1 - interval) / 2
  boundary <- as.numeric(x$reference$boundary_radius)[1L]
  rows <- lapply(groups, function(index) {
    dat <- tab[index, , drop = FALSE]
    values <- dat$niche_boundary_exceedance
    finite <- is.finite(values)
    q <- if (any(finite)) {
      stats::quantile(
        values[finite],
        probs = c(alpha, 0.5, 1 - alpha),
        names = FALSE,
        type = 8
      )
    } else {
      rep(NA_real_, 3L)
    }
    data.frame(
      cell = dat$cell[1L],
      time = dat$time[1L],
      scenario = dat$scenario[1L],
      scope = scope,
      n_models = sum(finite),
      model_agreement = if (sum(finite) >= 2L) {
        mean(values[finite] > tolerance)
      } else {
        NA_real_
      },
      lower_niche_boundary_exceedance = q[1L],
      median_niche_boundary_exceedance = q[2L],
      upper_niche_boundary_exceedance = q[3L],
      median_relative_exceedance = if (boundary > 0) q[2L] / boundary else NA,
      interval = interval,
      boundary_exceedance_tolerance = tolerance,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("climniche_model_agreement", "data.frame")
  out
}

#' Summarise interval changes in range level niche boundary exceedance
#'
#' @param x A `climniche_series` object.
#' @param metric Range summary used to quantify increase.
#' @param scope,aggregation_weight,area_weight,boundary_exceedance_tolerance
#'   Arguments passed to [climniche_range_summary()].
#'
#' @return A data frame with one row per model and scenario.
#'
#' @details
#' For a selected range summary \eqn{X_t}, interval rates are
#' \deqn{g_k = \frac{X_{t_{k+1}} - X_{t_k}}{t_{k+1} - t_k}.}
#' The function reports the largest positive \eqn{g_k} and its interval. Rates
#' should only be compared among analyses using compatible temporal units and
#' projection resolution.
#' @export
climniche_change_rate <- function(
    x,
    metric = c("exposed_fraction",
               "range_wide_relative_exceedance",
               "mean_niche_boundary_exceedance"),
    scope = c("current", "all"),
    aggregation_weight = NULL,
    area_weight = FALSE,
    boundary_exceedance_tolerance = NULL) {
  if (!inherits(x, "climniche_series")) {
    stop("x must be a climniche_series object.", call. = FALSE)
  }
  metric <- match.arg(metric)
  scope <- match.arg(scope)
  summary <- climniche_range_summary(
    x,
    scope = scope,
    aggregation_weight = aggregation_weight,
    area_weight = area_weight,
    boundary_exceedance_tolerance = boundary_exceedance_tolerance
  )
  groups <- split(
    seq_len(nrow(summary)),
    .series_group_key(summary$model, summary$scenario)
  )
  rows <- lapply(groups, function(index) {
    dat <- summary[index, , drop = FALSE]
    ord <- order(as.numeric(dat$time))
    dat <- dat[ord, , drop = FALSE]
    time_numeric <- .time_numeric(dat$time)
    value <- dat[[metric]]
    if (length(value) > 1L) {
      interval <- diff(time_numeric)
      change <- diff(value)
      rates <- change / interval
      positive <- which(is.finite(rates) & rates > 0)
    } else {
      change <- rates <- numeric()
      positive <- integer()
    }
    if (length(positive)) {
      selected <- positive[which.max(rates[positive])]
      start <- dat$time[selected]
      end <- dat$time[selected + 1L]
      maximum_change <- change[selected]
      maximum_rate <- rates[selected]
    } else {
      start <- dat$time[NA_integer_]
      end <- dat$time[NA_integer_]
      maximum_change <- 0
      maximum_rate <- 0
    }
    data.frame(
      model = dat$model[1L],
      scenario = dat$scenario[1L],
      scope = scope,
      metric = metric,
      n_projections = nrow(dat),
      first_time = dat$time[1L],
      last_time = dat$time[nrow(dat)],
      total_change = utils::tail(value, 1L) - value[1L],
      maximum_interval_increase = maximum_change,
      maximum_increase_rate = maximum_rate,
      interval_start = start,
      interval_end = end,
      time_unit = .time_unit(dat$time),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("climniche_change_rate", "data.frame")
  out
}

#' @export
summary.climniche_series <- function(
    object,
    scope = c("current", "all"),
    aggregation_weight = NULL,
    area_weight = FALSE,
    boundary_exceedance_tolerance = NULL,
    ...) {
  scope <- match.arg(scope)
  out <- list(
    projections = object$index,
    range = climniche_range_summary(
      object,
      scope = scope,
      aggregation_weight = aggregation_weight,
      area_weight = area_weight,
      boundary_exceedance_tolerance =
        boundary_exceedance_tolerance
    ),
    change_rate = climniche_change_rate(
      object,
      scope = scope,
      aggregation_weight = aggregation_weight,
      area_weight = area_weight,
      boundary_exceedance_tolerance =
        boundary_exceedance_tolerance
    )
  )
  class(out) <- "summary.climniche_series"
  out
}

#' @export
print.summary.climniche_series <- function(x, ...) {
  cat("Niche relative climate exposure through time\n\n")
  print(.series_display_table(x$range), row.names = FALSE)
  cat("\nMaximum interval increase\n")
  print(.series_display_table(x$change_rate), row.names = FALSE)
  invisible(x)
}
