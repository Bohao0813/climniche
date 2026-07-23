if (!exists("fit_climniche", mode = "function") &&
    requireNamespace("climniche", quietly = TRUE)) {
  library(climniche)
}

dynamic_internal <- function(name) {
  if (exists(name, mode = "function")) {
    return(get(name, mode = "function"))
  }
  getFromNamespace(name, "climniche")
}

dynamic_definitions <- dynamic_internal(".dynamic_definitions")()
stopifnot(identical(
  dynamic_definitions$name[seq_len(3)],
  c(
    "Weighted Niche Boundary Exceedance Fraction",
    "Conditional Relative Niche Boundary Exceedance",
    "Range Mean Relative Niche Boundary Exceedance"
  )
))
stopifnot(identical(
  dynamic_definitions$name[4:6],
  c(
    "Persistent Niche Boundary Exceedance Onset",
    "Time Weighted Niche Boundary Exceedance Fraction",
    "Cumulative Relative Niche Boundary Exceedance"
  )
))
stopifnot(all(grepl(
  "tolerance",
  dynamic_definitions$definition[c(1:5, 7)],
  ignore.case = TRUE
)))
stopifnot(!any(grepl(
  "departure|habitat",
  dynamic_definitions$name,
  ignore.case = TRUE
)))

set.seed(81)
dynamic_current <- matrix(rnorm(720), ncol = 6)
colnames(dynamic_current) <- paste0("climate_", seq_len(6))
rownames(dynamic_current) <- paste0("cell_", seq_len(nrow(dynamic_current)))
dynamic_weights <- seq(0.1, 1, length.out = nrow(dynamic_current))
dynamic_shift <- matrix(
  rep(c(0.35, -0.15, 0.20, 0.10, -0.08, 0.18),
      nrow(dynamic_current)),
  ncol = 6,
  byrow = TRUE
)
dynamic_future <- dynamic_current + dynamic_shift

# Projection against a stored reference must reproduce a direct fit.
dynamic_reference <- fit_climniche_reference(
  dynamic_current,
  occupied = dynamic_weights,
  sensitivity = seq(0.7, 1.3, length.out = 6),
  preprocess = FALSE
)
projected_fit <- project_climniche(dynamic_reference, dynamic_future)
stopifnot(isTRUE(all.equal(
  projected_fit$occupied_weight,
  dynamic_weights,
  tolerance = 0
)))
new_current <- dynamic_current
new_future <- dynamic_future
rownames(new_current) <- paste0("new_cell_", seq_len(nrow(new_current)))
rownames(new_future) <- rownames(new_current)
new_row_fit <- project_climniche(
  dynamic_reference,
  new_future,
  current = new_current
)
stopifnot(identical(new_row_fit$occupied_weight, rep(1, nrow(dynamic_current))))
direct_fit <- fit_climniche(
  dynamic_current,
  dynamic_future,
  occupied = dynamic_weights,
  sensitivity = seq(0.7, 1.3, length.out = 6),
  preprocess = FALSE
)
for (metric in c(
  "climate_change_amount",
  "niche_distance_change",
  "climate_reconfiguration",
  "niche_boundary_exceedance"
)) {
  stopifnot(isTRUE(all.equal(
    projected_fit[[metric]],
    direct_fit[[metric]],
    tolerance = 1e-12
  )))
}
stopifnot(isTRUE(all.equal(projected_fit$center, direct_fit$center,
                           tolerance = 0)))
stopifnot(isTRUE(all.equal(projected_fit$A, direct_fit$A, tolerance = 0)))

# Climatic Reconfiguration is the distance scaled angular term.
z0 <- sweep(projected_fit$current, 2L, projected_fit$center, "-")
z1 <- sweep(projected_fit$future, 2L, projected_fit$center, "-")
inner_product <- rowSums((z0 %*% projected_fit$A) * z1)
r0 <- projected_fit$niche_radius_current
r1 <- projected_fit$niche_radius_future
valid_angle <- r0 * r1 > sqrt(.Machine$double.eps)
cosine <- inner_product[valid_angle] / (r0[valid_angle] * r1[valid_angle])
angular_term <- 2 * r0[valid_angle] * r1[valid_angle] * (1 - cosine)
stopifnot(isTRUE(all.equal(
  projected_fit$climate_reconfiguration[valid_angle]^2,
  angular_term,
  tolerance = 1e-10,
  check.attributes = FALSE
)))

fractions <- c(0.2, 0.5, 0.8, 1)
dynamic_projections <- lapply(fractions, function(fraction) {
  dynamic_current + fraction * dynamic_shift
})
dynamic_series <- fit_climniche_series(
  current = dynamic_current,
  future = dynamic_projections,
  time = c(2030, 2050, 2070, 2090),
  occupied = dynamic_weights,
  sensitivity = seq(0.7, 1.3, length.out = 6),
  preprocess = FALSE
)
stopifnot(inherits(dynamic_series, "climniche_series"))
stopifnot(all(vapply(dynamic_series$fits, function(fit) {
  isTRUE(all.equal(fit$center, dynamic_series$reference$center, tolerance = 0))
}, logical(1))))
stopifnot(all(vapply(dynamic_series$fits, function(fit) {
  isTRUE(all.equal(fit$boundary_radius,
                   dynamic_series$reference$boundary_radius,
                   tolerance = 0))
}, logical(1))))
pooled_distance_change <- unlist(lapply(
  dynamic_series$fits,
  function(fit) abs(fit$niche_distance_change)
), use.names = FALSE)
expected_tolerance <- as.numeric(stats::quantile(
  pooled_distance_change,
  probs = 0.10,
  names = FALSE,
  type = 8
))
stopifnot(isTRUE(all.equal(
  dynamic_series$descriptor_settings$tolerance,
  expected_tolerance,
  tolerance = 0
)))
stopifnot(all(vapply(dynamic_series$fits, function(fit) {
  identical(fit$descriptor_settings, dynamic_series$descriptor_settings)
}, logical(1))))

# The range mean equals exposed fraction multiplied by conditional severity.
range_summary <- climniche_range_summary(dynamic_series)
stopifnot(isTRUE(all.equal(
  range_summary$range_wide_relative_exceedance,
  range_summary$exposed_fraction *
    range_summary$conditional_relative_exceedance,
  tolerance = 1e-12
)))

# The Mediterranean time-series tables retain the same range decomposition.
case_path <- system.file(
  "extdata",
  "mediterranean_anchovy",
  package = "climniche"
)
case_range <- utils::read.csv(file.path(
  case_path,
  "anchovy_climniche_time_range_summary.csv"
))
case_departure <- utils::read.csv(file.path(
  case_path,
  "anchovy_climniche_time_departure_summary.csv"
))
stopifnot(identical(
  as.integer(case_range$time),
  c(2030L, 2050L, 2070L, 2090L)
))
stopifnot(isTRUE(all.equal(
  case_range$range_wide_relative_exceedance,
  case_range$exposed_fraction *
    case_range$conditional_relative_exceedance,
  tolerance = 1e-12
)))
stopifnot(
  nrow(case_departure) == 1L,
  case_departure$proportion_with_persistent_departure >= 0,
  case_departure$proportion_with_persistent_departure <= 1
)
weighted_range <- climniche_range_summary(
  dynamic_series,
  aggregation_weight = seq_len(nrow(dynamic_current))
)
stopifnot(!isTRUE(all.equal(
  range_summary$aggregation_weight_sum,
  weighted_range$aggregation_weight_sum
)))
zero_range <- climniche_range_summary(
  dynamic_series,
  boundary_exceedance_tolerance = 1e6
)
stopifnot(all(zero_range$exposed_fraction == 0))
stopifnot(all(zero_range$conditional_relative_exceedance == 0))
stopifnot(all(zero_range$range_wide_relative_exceedance == 0))

# Persistent exceedance is tied to projection times, and a subsequent return
# below the boundary is counted when a qualifying run ends.
departure_series <- dynamic_series
departure_values <- c(0, 0.2, 0.3, 0)
for (i in seq_along(departure_series$fits)) {
  departure_series$fits[[i]]$niche_boundary_exceedance[1] <-
    departure_values[i]
  departure_series$fits[[i]]$outside_niche_exceedance[1] <-
    departure_values[i]
}
departure <- climniche_departure(
  departure_series,
  scope = "all",
  persistence = 2,
  boundary_exceedance_tolerance = 0
)
departure_cell <- departure[departure$cell == 1, , drop = FALSE]
stopifnot(identical(as.numeric(departure_cell$first_persistent_departure),
                    2050))
stopifnot(identical(departure_cell$reentry_count, 1L))
stopifnot(isTRUE(departure_cell$reentered))

date_series <- departure_series
date_series$index$time <- as.Date(c(
  "2030-01-01", "2050-01-01", "2070-01-01", "2090-01-01"
))
date_departure <- climniche_departure(
  date_series,
  scope = "all",
  persistence = 2,
  boundary_exceedance_tolerance = 0
)
stopifnot(inherits(date_departure$first_persistent_departure, "Date"))
stopifnot(identical(unique(date_departure$time_unit), "days"))

posix_series <- departure_series
posix_series$index$time <- as.POSIXct(
  c("2030-01-01", "2050-01-01", "2070-01-01", "2090-01-01"),
  tz = "UTC"
)
posix_departure <- climniche_departure(
  posix_series,
  scope = "all",
  persistence = 2,
  boundary_exceedance_tolerance = 0
)
stopifnot(identical(unique(posix_departure$time_unit), "seconds"))

# Agreement is the fraction of supplied models beyond the boundary.
ensemble_series <- fit_climniche_series(
  current = dynamic_current,
  future = list(dynamic_future, dynamic_future),
  time = c(2050, 2050),
  model = c("model_a", "model_b"),
  scenario = "SSP2-4.5",
  occupied = dynamic_weights,
  sensitivity = seq(0.7, 1.3, length.out = 6),
  preprocess = FALSE
)
ensemble_series$fits[[1]]$niche_boundary_exceedance[1] <- 0.2
ensemble_series$fits[[1]]$outside_niche_exceedance[1] <- 0.2
ensemble_series$fits[[2]]$niche_boundary_exceedance[1] <- 0
ensemble_series$fits[[2]]$outside_niche_exceedance[1] <- 0
agreement <- climniche_model_agreement(
  ensemble_series,
  scope = "all",
  boundary_exceedance_tolerance = 0
)
stopifnot(isTRUE(all.equal(
  agreement$model_agreement[agreement$cell == 1],
  0.5
)))
single_model_agreement <- climniche_model_agreement(dynamic_series)
stopifnot(all(is.na(single_model_agreement$model_agreement)))
ensemble_report <- climniche_series_report(ensemble_series, scope = "all")
stopifnot(!is.null(ensemble_report$model_agreement_summary))

change_rate <- climniche_change_rate(dynamic_series)
stopifnot(nrow(change_rate) == 1L)
stopifnot(is.finite(change_rate$maximum_increase_rate))

dynamic_report <- climniche_series_report(dynamic_series)
stopifnot(inherits(dynamic_report, "climniche_series_report"))
stopifnot(is.null(dynamic_report$model_agreement_summary))
date_report <- climniche_series_report(
  date_series,
  scope = "all",
  persistence = 2,
  boundary_exceedance_tolerance = 0
)
stopifnot(inherits(
  date_report$departure_summary$median_first_persistent_departure,
  "Date"
))
report_file <- tempfile(fileext = ".md")
write_climniche_series_report(dynamic_report, report_file)
stopifnot(file.exists(report_file), file.info(report_file)$size > 0)
unlink(report_file)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  time_plot <- plot_climniche_time(dynamic_series)
  stopifnot(inherits(time_plot, "ggplot"))
  stopifnot(is.null(time_plot$scales$get_scales("fill")))
  stopifnot(identical(
    time_plot$labels$y,
    "Range Mean Relative Niche Boundary Exceedance"
  ))
  fraction_plot <- plot_climniche_time(
    dynamic_series,
    metric = "exposed_fraction"
  )
  stopifnot(identical(
    fraction_plot$labels$y,
    "Weighted Niche Boundary Exceedance Fraction"
  ))
}

if (requireNamespace("raster", quietly = TRUE)) {
  r1 <- raster::raster(nrows = 3, ncols = 4, xmn = -2, xmx = 2,
                       ymn = 40, ymx = 43)
  r2 <- r1
  raster::values(r1) <- seq(-1, 1, length.out = raster::ncell(r1))
  raster::values(r2) <- rep(c(-0.8, -0.2, 0.4), length.out = raster::ncell(r2))
  raster_current <- raster::stack(r1, r2)
  names(raster_current) <- c("temperature", "oxygen")
  raster_future <- lapply(c(0.2, 0.5, 0.8), function(fraction) {
    raster_current + fraction * 0.3
  })
  missing_values <- raster::getValues(raster_future[[2]])
  missing_values[1, 1] <- NA_real_
  raster::values(raster_future[[2]]) <- missing_values
  raster_weights <- r1
  expected_weights <- seq(0.1, 1, length.out = raster::ncell(r1))
  raster::values(raster_weights) <- expected_weights
  raster_series <- fit_climniche_series(
    current = raster_current,
    future = raster_future,
    time = c(2030, 2060, 2090),
    occupied = raster_weights,
    occupied_threshold = 0.3,
    preprocess = FALSE
  )
  full_expected_weights <- expected_weights
  full_expected_weights[full_expected_weights <= 0.3] <- 0
  expected_weights <- full_expected_weights[-1]
  stopifnot(isTRUE(all.equal(
    raster_series$fits[[1]]$occupied_weight,
    expected_weights,
    tolerance = 0
  )))
  stopifnot(isTRUE(all.equal(
    raster_series$reference$reference_weight,
    full_expected_weights,
    tolerance = 0
  )))
  raster_reference_fit <- fit_climniche_raster(
    current = raster_current,
    future = raster_current,
    occupied = raster_weights,
    occupied_threshold = 0.3,
    preprocess = FALSE,
    tolerance = 0
  )
  stopifnot(isTRUE(all.equal(
    raster_series$reference$center,
    raster_reference_fit$center,
    tolerance = 0
  )))
  stopifnot(isTRUE(all.equal(
    raster_series$reference$boundary_radius,
    raster_reference_fit$boundary_radius,
    tolerance = 0
  )))
  stopifnot(length(unique(vapply(
    raster_series$fits,
    function(fit) paste(fit$cell_index, collapse = ","),
    character(1)
  ))) == 1L)
  stopifnot(identical(raster_series$fits[[1]]$cell_index, 2:12))
  area_summary <- climniche_range_summary(
    raster_series,
    area_weight = TRUE
  )
  stopifnot(all(is.finite(area_summary$aggregation_weight_sum)))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    departure_map <- plot_climniche_departure_map(
      raster_series,
      metric = "departure_time_fraction",
      scope = "all",
      title = FALSE,
      legend_title = FALSE
    )
    stopifnot(inherits(departure_map, "ggplot"))
    raster_ensemble <- fit_climniche_series(
      current = raster_current,
      future = list(raster_future[[1]], raster_future[[3]]),
      time = c(2050, 2050),
      model = c("model_a", "model_b"),
      occupied = raster_weights,
      occupied_threshold = 0.3,
      preprocess = FALSE
    )
    agreement_map <- plot_climniche_departure_map(
      raster_ensemble,
      metric = "model_agreement",
      time = 2050,
      scope = "all",
      title = FALSE,
      legend_title = FALSE
    )
    stopifnot(inherits(agreement_map, "ggplot"))
  }
}

if (requireNamespace("terra", quietly = TRUE)) {
  terra_current <- terra::rast(nrows = 3, ncols = 4, nlyrs = 2,
                               xmin = -2, xmax = 2, ymin = 40, ymax = 43)
  terra::values(terra_current) <- cbind(
    seq(-1, 1, length.out = terra::ncell(terra_current)),
    rep(c(-0.8, -0.2, 0.4), length.out = terra::ncell(terra_current))
  )
  names(terra_current) <- c("temperature", "oxygen")
  terra_future <- lapply(c(0.2, 0.5), function(fraction) {
    terra_current + fraction * 0.3
  })
  missing_values <- terra::values(terra_future[[2]], mat = TRUE)
  missing_values[nrow(missing_values), 2] <- NA_real_
  terra::values(terra_future[[2]]) <- missing_values
  terra_weights <- terra_current[[1]]
  expected_weights <- seq(0.1, 1, length.out = terra::ncell(terra_weights))
  terra::values(terra_weights) <- expected_weights
  terra_series <- fit_climniche_series(
    current = terra_current,
    future = terra_future,
    time = c(2050, 2090),
    occupied = terra_weights,
    occupied_threshold = 0.3,
    preprocess = FALSE
  )
  full_expected_weights <- expected_weights
  full_expected_weights[full_expected_weights <= 0.3] <- 0
  expected_weights <- full_expected_weights[-length(full_expected_weights)]
  stopifnot(isTRUE(all.equal(
    terra_series$fits[[1]]$occupied_weight,
    expected_weights,
    tolerance = 0
  )))
  stopifnot(isTRUE(all.equal(
    terra_series$reference$reference_weight,
    full_expected_weights,
    tolerance = 0
  )))
  terra_reference_fit <- fit_climniche_terra(
    current = terra_current,
    future = terra_current,
    occupied = terra_weights,
    occupied_threshold = 0.3,
    preprocess = FALSE,
    tolerance = 0
  )
  stopifnot(isTRUE(all.equal(
    terra_series$reference$center,
    terra_reference_fit$center,
    tolerance = 0
  )))
  stopifnot(isTRUE(all.equal(
    terra_series$reference$boundary_radius,
    terra_reference_fit$boundary_radius,
    tolerance = 0
  )))
}
