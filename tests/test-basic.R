if (!exists("fit_climniche", mode = "function") &&
    requireNamespace("climniche", quietly = TRUE)) {
  library(climniche)
}

sim <- simulate_climniche(n = 500, seed = 10)

toward <- fit_climniche(
  current = sim$current,
  future = sim$future_toward,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity,
  scale = FALSE
)

away <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity,
  scale = FALSE
)

away_numeric_occupied <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = as.numeric(sim$occupied),
  occupied_threshold = 0.5,
  center = sim$center,
  sensitivity = sim$sensitivity,
  scale = FALSE
)

stopifnot(inherits(toward, "climniche_fit"))
stopifnot(identical(away$occupied, away_numeric_occupied$occupied))
stopifnot(length(toward$climate_change_amount) == nrow(sim$current))
stopifnot(length(toward$climate_reconfiguration) == nrow(sim$current))
stopifnot(length(toward$niche_boundary_exceedance) == nrow(sim$current))
stopifnot(is.null(away$classification))
stopifnot(is.null(away$mixed_variable_response))
stopifnot(is.null(away$classification_settings))
stopifnot(isTRUE(all.equal(
  toward$climate_reconfiguration,
  toward$composition_change,
  tolerance = 0
)))
stopifnot(isTRUE(all.equal(
  toward$niche_boundary_exceedance,
  toward$outside_niche_exceedance,
  tolerance = 0
)))
fit_summary <- summary(toward)
stopifnot(isTRUE(all.equal(
  fit_summary$climate_reconfiguration,
  fit_summary$composition_change
)))
stopifnot(isTRUE(all.equal(
  fit_summary$niche_boundary_exceedance,
  fit_summary$outside_niche_exceedance
)))
stopifnot(isTRUE(all.equal(
  fit_summary$boundary_distance,
  toward$boundary_distance
)))
summary_text <- capture.output(print(fit_summary))
stopifnot(any(grepl("Climatic Reconfiguration", summary_text, fixed = TRUE)))
stopifnot(!any(grepl("composition_change", summary_text, fixed = TRUE)))
stopifnot(all(toward$climate_reconfiguration >= -sqrt(.Machine$double.eps)))
stopifnot(all(abs(toward$change_alignment[!is.na(toward$change_alignment)]) <=
                1 + sqrt(.Machine$double.eps)))
stopifnot(isTRUE(all.equal(
  toward$climate_change_amount^2,
  toward$niche_distance_change^2 + toward$climate_reconfiguration^2,
  tolerance = 1e-7
)))
stopifnot(mean(toward$niche_distance_change) < mean(away$niche_distance_change))
stopifnot(sum(away$niche_boundary_exceedance > 0) >=
             sum(toward$niche_boundary_exceedance > 0))

# Controlled paths distinguish radial, directional and boundary terms.
geometry_current <- matrix(
  rep(c(1, 0), 5),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(
    c("unchanged", "toward", "away", "around", "cross_centre"),
    c("climate_1", "climate_2")
  )
)
geometry_future <- rbind(
  unchanged = c(1, 0),
  toward = c(0, 0),
  away = c(2, 0),
  around = c(0, 1),
  cross_centre = c(-1, 0)
)
colnames(geometry_future) <- colnames(geometry_current)
geometry_fit <- fit_climniche(
  geometry_current,
  geometry_future,
  occupied = rep(1, nrow(geometry_current)),
  center = c(0, 0),
  A = diag(2),
  boundary = 0.95,
  scale = FALSE,
  preprocess = FALSE,
  tolerance = 0
)
stopifnot(isTRUE(all.equal(
  unname(geometry_fit$climate_change_amount),
  c(0, 1, 1, sqrt(2), 2),
  tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  unname(geometry_fit$niche_distance_change),
  c(0, -1, 1, 0, 0),
  tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  unname(geometry_fit$climate_reconfiguration),
  c(0, 0, 0, sqrt(2), 2),
  tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  unname(geometry_fit$niche_boundary_exceedance),
  c(0, 0, 1, 0, 0),
  tolerance = 1e-12
)))

# Exceedance is a future radial state, not a crossing indicator.
baseline_fit <- fit_climniche(
  sim$current,
  sim$current,
  occupied = sim$occupied,
  sensitivity = sim$sensitivity,
  boundary = 0.95
)
stopifnot(all(baseline_fit$climate_change_amount == 0))
stopifnot(any(
  baseline_fit$niche_boundary_exceedance[baseline_fit$occupied] > 0
))

# Weighted percentiles include the upper endpoint and tied reference values.
percentile_endpoints <- niche_percentile(
  psi_current = c(1, 2, 3),
  psi_future = c(1, 2, 3),
  occupied = c(1, 1, 1)
)
stopifnot(isTRUE(all.equal(
  percentile_endpoints$current,
  c(1 / 3, 2 / 3, 1),
  tolerance = 1e-12
)))
percentile_ties <- niche_percentile(
  psi_current = c(1, 2, 2, 3),
  psi_future = c(2, 3, 1, 2),
  occupied = c(1, 2, 3, 4)
)
stopifnot(isTRUE(all.equal(
  percentile_ties$future,
  c(0.6, 1, 0.1, 0.6),
  tolerance = 1e-12
)))
single_percentile <- niche_percentile(1, 1, occupied = 2)
stopifnot(identical(single_percentile$current, 1))
stopifnot(identical(single_percentile$future, 1))

# Numerical negative eigenvalues within tolerance are projected to zero.
near_psd_fit <- fit_climniche(
  current = matrix(c(0, 0), nrow = 1),
  future = matrix(c(0, 1), nrow = 1),
  occupied = 1,
  center = c(0, 0),
  A = diag(c(1, -1e-9)),
  scale = FALSE,
  preprocess = FALSE
)
stopifnot(min(eigen(near_psd_fit$A, symmetric = TRUE)$values) >= -1e-14)
stopifnot(isTRUE(all.equal(
  near_psd_fit$climate_change_amount^2,
  near_psd_fit$niche_distance_change^2 +
    near_psd_fit$climate_reconfiguration^2,
  tolerance = 1e-12
)))

tab <- climniche_table(away)
stopifnot(all(c("occupied_weight", "climate_reconfiguration",
                "niche_boundary_exceedance", "composition_change",
                "outside_niche_exceedance", "radial_direction",
                "boundary_status") %in% names(tab)))
stopifnot(!("class" %in% names(tab)))
stopifnot(!("mixed_variable_response" %in% names(tab)))
stopifnot(identical(levels(away$radial_direction), c(
  "Toward realised niche centre",
  "Limited Niche Distance Shift",
  "Away from realised niche centre"
)))
stopifnot(identical(levels(away$boundary_status), c(
  "Within empirical niche boundary",
  "Beyond empirical niche boundary"
)))
stopifnot(is.finite(away$descriptor_settings$tolerance))
stopifnot(isTRUE(all.equal(
  away$descriptor_settings$tolerance,
  away$threshold_settings$tolerance
)))

custom_descriptor <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity,
  scale = FALSE,
  tolerance = 1e6,
  boundary_exceedance_tolerance = 1e6
)
stopifnot(all(custom_descriptor$radial_direction ==
                "Limited Niche Distance Shift"))
stopifnot(all(custom_descriptor$boundary_status ==
                "Within empirical niche boundary"))
stopifnot(isTRUE(all.equal(custom_descriptor$descriptor_settings$tolerance,
                           1e6)))
stopifnot(isTRUE(all.equal(
  custom_descriptor$descriptor_settings$boundary_exceedance_tolerance,
  1e6
)))

simple_current <- matrix(c(
  0, 0,
  10, 0,
  0, 10,
  10, 10
), ncol = 2, byrow = TRUE)
simple_future <- simple_current + 0.2
simple_weights <- c(1, 10, 0, 0)
weighted_fit <- fit_climniche(
  simple_current,
  simple_future,
  occupied = simple_weights,
  sensitivity = c(1, 1),
  scale = FALSE
)
expected_center <- colSums(simple_current * simple_weights) /
  sum(simple_weights)
stopifnot(isTRUE(all.equal(unname(weighted_fit$center), expected_center)))
stopifnot(!isTRUE(all.equal(
  unname(weighted_fit$center),
  colMeans(simple_current[1:2, ])
)))

threshold_fit <- fit_climniche(
  simple_current,
  simple_future,
  occupied = c(0.2, 0.8, 2, 0),
  occupied_threshold = 0.5,
  sensitivity = c(1, 1),
  scale = FALSE
)
stopifnot(isTRUE(all.equal(threshold_fit$occupied_weight,
                           c(0, 0.8, 2, 0))))
stopifnot(isTRUE(all.equal(threshold_fit$occupied_weight[threshold_fit$occupied],
                           c(0.8, 2))))

index_fit <- fit_climniche(
  simple_current,
  simple_future,
  occupied = c(1, 3),
  sensitivity = c(1, 1),
  scale = FALSE
)
stopifnot(isTRUE(all.equal(index_fit$occupied_weight, c(1, 0, 1, 0))))

length_n_weight_fit <- fit_climniche(
  simple_current,
  simple_future,
  occupied = c(1, 2, 3, 4),
  sensitivity = c(1, 1),
  scale = FALSE
)
expected_length_n_center <- colSums(simple_current * c(1, 2, 3, 4)) / 10
stopifnot(isTRUE(all.equal(unname(length_n_weight_fit$center),
                            expected_length_n_center)))

summary_future <- simple_current
summary_future[, 1] <- summary_future[, 1] + c(0, 1, 2, 3)
summary_weight <- c(100, 1, 0, 0)
summary_fit <- fit_climniche(
  simple_current,
  summary_future,
  occupied = summary_weight,
  sensitivity = c(1, 1),
  scale = FALSE,
  preprocess = FALSE
)
summary_current <- summary(summary_fit)
summary_all <- summary(summary_fit, scope = "all")
expected_weighted_amount <- sum(
  summary_fit$climate_change_amount * summary_weight
) / sum(summary_weight)
stopifnot(identical(summary_current$scope, "current"))
stopifnot(identical(summary_all$scope, "all"))
stopifnot(isTRUE(all.equal(
  as.numeric(summary_current$climate_change_amount["Mean"]),
  expected_weighted_amount,
  tolerance = 1e-12
)))
stopifnot(!isTRUE(all.equal(
  as.numeric(summary_current$climate_change_amount["Mean"]),
  as.numeric(summary_all$climate_change_amount["Mean"])
)))
plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
plotted_current <- plot(summary_fit, type = "amount")
plotted_all <- plot(summary_fit, type = "reconfiguration", scope = "all")
grDevices::dev.off()
unlink(plot_file)
stopifnot(inherits(plotted_current, "climniche_fit"))
stopifnot(inherits(plotted_all, "climniche_fit"))

negative_weight <- try(fit_climniche(
  simple_current,
  simple_future,
  occupied = c(1, -1, 0, 0),
  sensitivity = c(1, 1),
  scale = FALSE
), silent = TRUE)
stopifnot(inherits(negative_weight, "try-error"))

pre_current <- cbind(
  mean_temperature = c(0, 1, 2, 3, 4, 5),
  duplicate_temperature = c(0, 2, 4, 6, 8, 10),
  salinity_range = c(2, 1, 2, 1, 2, 1),
  constant_layer = rep(7, 6)
)
pre_future <- pre_current + cbind(
  mean_temperature = rep(0.2, 6),
  duplicate_temperature = rep(0.4, 6),
  salinity_range = rep(0.1, 6),
  constant_layer = rep(0, 6)
)
pre_fit <- fit_climniche(
  pre_current,
  pre_future,
  occupied = rep(1, 6),
  sensitivity = c(1, 1, 1, 1),
  preprocess_correlation = 0.95
)
stopifnot(isTRUE(pre_fit$preprocessing$settings$enabled))
stopifnot(sum(pre_fit$preprocessing$retained_variables %in%
                c("mean_temperature", "duplicate_temperature")) == 1L)
stopifnot("salinity_range" %in% pre_fit$preprocessing$retained_variables)
stopifnot("constant_layer" %in%
            pre_fit$preprocessing$removed_variables$variable)
stopifnot(ncol(pre_fit$current) == 2L)
stopifnot(is.null(pre_fit$preprocessing$current))
stopifnot(is.null(pre_fit$preprocessing$future))
stopifnot(isTRUE(all.equal(
  pre_fit$climate_change_amount^2,
  pre_fit$niche_distance_change^2 + pre_fit$climate_reconfiguration^2,
  tolerance = 1e-7
)))
pre_fit_rescaled <- fit_climniche(
  pre_current * 1000,
  pre_future * 1000,
  occupied = rep(1, 6),
  sensitivity = c(1, 1, 1, 1),
  preprocess_correlation = 0.95
)
stopifnot(identical(
  pre_fit$preprocessing$retained_variables,
  pre_fit_rescaled$preprocessing$retained_variables
))

pre_off <- fit_climniche(
  pre_current,
  pre_future,
  occupied = rep(1, 6),
  sensitivity = c(1, 1, 1, 1),
  scale = FALSE,
  preprocess = FALSE
)
stopifnot(!isTRUE(pre_off$preprocessing$settings$enabled))
stopifnot(ncol(pre_off$current) == 4L)

pre_A <- fit_climniche(
  pre_current,
  pre_future,
  occupied = rep(1, 6),
  A = diag(c(1, 2, 3, 4)),
  preprocess_correlation = 0.95
)
stopifnot(identical(dim(pre_A$A), c(2L, 2L)))

cnfa_like <- list(
  mf = rep(0, 4),
  sf = c(1, 1, 1, 1),
  co = diag(4),
  eig = rep(1, 4)
)
pre_cnfa <- fit_climniche(
  pre_current,
  pre_future,
  occupied = rep(1, 6),
  cnfa = cnfa_like,
  metric = "factor",
  scale = FALSE,
  preprocess = FALSE
)
stopifnot(inherits(pre_cnfa, "climniche_fit"))
stopifnot(ncol(pre_cnfa$current) == 4L)

bad_cnfa_preprocess <- try(fit_climniche(
  pre_current,
  pre_future,
  occupied = rep(1, 6),
  cnfa = cnfa_like,
  metric = "factor",
  scale = FALSE
), silent = TRUE)
stopifnot(inherits(bad_cnfa_preprocess, "try-error"))

bad_cnfa_scale <- try(fit_climniche(
  pre_current,
  pre_future,
  occupied = rep(1, 6),
  cnfa = cnfa_like,
  metric = "factor",
  preprocess = FALSE
), silent = TRUE)
stopifnot(inherits(bad_cnfa_scale, "try-error"))

A_full <- diag(seq(1, 1.5, length.out = ncol(sim$current)))
A_full[1, 2] <- A_full[2, 1] <- 0.25
fit_full <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  A = A_full,
  scale = FALSE
)
stopifnot(isTRUE(all.equal(rowSums(fit_full$variable_contribution),
                           fit_full$psi_future - fit_full$psi_current,
                           tolerance = 1e-7)))
bad_A <- A_full
bad_A[1, 2] <- 2
bad_A[2, 1] <- 0
bad_fit <- try(fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  A = bad_A,
  scale = FALSE
), silent = TRUE)
stopifnot(inherits(bad_fit, "try-error"))

summary_data <- climniche_summary_figure_data(away, max_points = 80)
stopifnot(inherits(summary_data, "climniche_summary_figure_data"))
stopifnot(nrow(summary_data$plane) <= 80)
stopifnot(all(c("x_mid", "y_mid", "total_weight") %in%
                names(summary_data$plane_bins)))
stopifnot(all(c("x_mid", "y_mid", "total_weight") %in%
                names(summary_data$reconfiguration_bins)))
stopifnot(all(c("metric", "value", "weight") %in% names(summary_data$metrics)))
stopifnot(!("classes" %in% names(summary_data)))
stopifnot(all(c("descriptor", "level", "proportion") %in%
                names(summary_data$descriptors)))
descriptor_counts <- table(summary_data$descriptors$descriptor)
stopifnot(identical(
  as.integer(descriptor_counts["Niche Distance Shift direction"]),
  3L
))
stopifnot(identical(
  as.integer(descriptor_counts["Niche boundary status"]),
  2L
))
stopifnot(all(c("metric", "x", "width", "proportion") %in%
                names(summary_data$metric_histograms)))
stopifnot("mean_absolute_share" %in% names(summary_data$variables))
stopifnot(all(c("boundary_quantile", "prop_exceeded") %in%
                names(summary_data$boundary)))

report <- climniche_report(away, scope = "current")
formal_metric_names <- c(
  "Climatic Displacement",
  "Niche Distance Shift",
  "Climatic Reconfiguration",
  "Niche Boundary Exceedance"
)
stopifnot(all(formal_metric_names %in% report$metric_definitions$metric))
stopifnot(all(formal_metric_names %in% report$metric_summary$metric))
stopifnot(nrow(report$metric_summary) == 12L)
stopifnot(all(table(report$metric_summary$metric) == 3L))
stopifnot(is.null(report$class_summary))
stopifnot(is.null(report$classification_settings))
stopifnot(nrow(report$descriptor_summary) == 5L)
tmp_report <- tempfile(fileext = ".md")
write_climniche_report(report, tmp_report)
report_text <- paste(readLines(tmp_report, warn = FALSE), collapse = "\n")
stopifnot(all(vapply(formal_metric_names, grepl, logical(1),
                     x = report_text, fixed = TRUE)))

summ <- climniche_summary(away)
stopifnot(!("prop_stable" %in% names(summ)))
stopifnot(!("prop_niche_exceedance" %in% names(summ)))
stopifnot(all(c("q90_niche_distance_change",
                "q90_climate_reconfiguration",
                "boundary_exceedance_tolerance") %in% names(summ)))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_new <- plot_climniche_distribution(away,
                                       metric = "climate_reconfiguration")
  p_old <- plot_climniche_distribution(away,
                                       metric = "composition_change")
  stopifnot(inherits(p_new, "ggplot"))
  stopifnot(inherits(p_old, "ggplot"))
  stopifnot(identical(p_new$labels$title, "Climatic Reconfiguration"))
  stopifnot(identical(p_new$labels$x, "Climatic Reconfiguration"))
  stopifnot(isTRUE(all.equal(p_new$data$value, p_old$data$value)))
  p_exposure <- plot_climniche_exposure(away)
  bad_colour <- try(plot_climniche_exposure(
    away, colour_by = "classification"
  ), silent = TRUE)
  stopifnot(inherits(p_exposure, "ggplot"))
  stopifnot(inherits(bad_colour, "try-error"))
  stopifnot(identical(p_exposure$labels$x, "Climatic Displacement"))
  stopifnot(identical(p_exposure$labels$y, "Niche Distance Shift"))
  summary_figure_plot <- plot_climniche_summary_figure(away, max_points = 80)
  stopifnot(inherits(summary_figure_plot, "patchwork") ||
              is.list(summary_figure_plot))
}

if (requireNamespace("raster", quietly = TRUE)) {
  single_current <- raster::raster(
    nrows = 4, ncols = 4, xmn = 0, xmx = 4, ymn = 0, ymx = 4,
    crs = "+proj=longlat +datum=WGS84"
  )
  single_values <- seq_len(raster::ncell(single_current))
  single_values[1] <- Inf
  raster::values(single_current) <- single_values
  single_future <- single_current
  future_values <- seq_len(raster::ncell(single_future)) + 0.5
  raster::values(single_future) <- future_values
  single_fit <- fit_climniche_raster(single_current, single_future)
  stopifnot(ncol(single_fit$current) == 1L)
  stopifnot(nrow(single_fit$current) == 15L)
  stopifnot(!single_fit$raster_complete[1])

  r1 <- raster::raster(nrows = 4, ncols = 4, xmn = 0, xmx = 4,
                       ymn = 0, ymx = 4, crs = "+proj=longlat +datum=WGS84")
  raster::values(r1) <- seq_len(raster::ncell(r1))
  r2 <- r1
  raster::values(r2) <- rev(seq_len(raster::ncell(r2)))
  cur <- raster::stack(r1, r2)
  fut <- raster::stack(r1 + 0.5, r2 + 1)
  names(cur) <- c("temperature", "salinity")
  names(fut) <- c("temperature", "salinity")
  occ <- r1
  raster::values(occ) <- c(rep(0, 7), 0.5, rep(0.7, 4), rep(1, 4))
  domain <- r1
  raster::values(domain) <- c(rep(0, 4), rep(1, 12))

  rf <- fit_climniche_raster(cur, fut, occupied = occ,
                             occupied_threshold = 0.5)
  stopifnot(inherits(rf, "climniche_fit"))
  stopifnot(length(rf$occupied) == 8)
  stopifnot(!any(rf$occupied_weight == 0.5))
  stopifnot(isTRUE(all.equal(sort(unique(rf$occupied_weight[rf$occupied])),
                             c(0.7, 1))))
  stopifnot(methods::is(rf$rasters$climate_change_amount, "RasterLayer"))
  stopifnot(!is.null(rf$rasters$climate_reconfiguration))
  stopifnot(!is.null(rf$rasters$niche_boundary_exceedance))
  stopifnot(identical(names(rf$rasters$climate_change_amount),
                      "climate_change_amount"))
  stopifnot(is.null(rf$rasters$classification))
  stopifnot(is.null(rf$class_lookup))
  stopifnot(identical(rf$call[[1]], as.name("fit_climniche_raster")))
  stopifnot(identical(climniche_table(rf, scope = "all")$cell,
                      which(rf$raster_complete)))
  rf_reordered <- fit_climniche_raster(
    cur, fut[[2:1]], occupied = occ, occupied_threshold = 0.5
  )
  stopifnot(isTRUE(all.equal(
    rf$climate_change_amount,
    rf_reordered$climate_change_amount,
    tolerance = 0
  )))
  bad_raster_occupied <- try(fit_climniche_raster(
    cur, fut, occupied = raster::stack(occ, occ)
  ), silent = TRUE)
  stopifnot(inherits(bad_raster_occupied, "try-error"))
  rf_domain <- fit_climniche_raster(cur, fut, occupied = occ,
                                    occupied_threshold = 0.5,
                                    domain = domain)
  stopifnot(length(rf_domain$climate_change_amount) == 12)
  stopifnot(length(rf_domain$occupied) == 8)
  stopifnot(isTRUE(all.equal(sort(unique(rf_domain$occupied_weight[rf_domain$occupied])),
                             c(0.7, 1))))

  # Future missing values do not alter the current reference centre.
  reference_current <- raster::raster(
    nrows = 2, ncols = 2, xmn = 0, xmx = 2, ymn = 0, ymx = 2
  )
  raster::values(reference_current) <- c(1, 2, 3, 100)
  reference_future <- reference_current + 1
  reference_future_values <- raster::getValues(reference_future)
  reference_future_values[4] <- NA_real_
  raster::values(reference_future) <- reference_future_values
  reference_weight <- reference_current
  raster::values(reference_weight) <- c(1, 0, 0, 100)
  reference_fit <- fit_climniche_raster(
    reference_current,
    reference_future,
    occupied = reference_weight,
    scale = FALSE,
    preprocess = FALSE
  )
  stopifnot(isTRUE(all.equal(
    unname(reference_fit$center),
    (1 + 100 * 100) / 101,
    tolerance = 1e-12
  )))
  stopifnot(identical(reference_fit$cell_index, 1:3))
  stopifnot(identical(reference_fit$reference_cell_index, 1:4))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_map_new <- plot_climniche_map(rf, metric = "niche_boundary_exceedance",
                                    occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_map_old <- plot_climniche_map(rf, metric = "outside_niche_exceedance",
                                    occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_maps <- plot_climniche_maps(
      rf, occupied = occ, occupied_only = TRUE,
      occupied_threshold = 0.5, degree_labels = "hemisphere"
    )
    stopifnot(inherits(p_map_new, "ggplot"))
    stopifnot(inherits(p_map_old, "ggplot"))
    stopifnot(!identical(
      p_map_new$scales$get_scales("fill")$guide,
      "none"
    ))
    p_map_hidden_legend <- plot_climniche_map(
      rf,
      metric = "niche_boundary_exceedance",
      show_legend = FALSE
    )
    stopifnot(identical(
      p_map_hidden_legend$theme$legend.position,
      "none"
    ))
    stopifnot(inherits(p_maps, "patchwork") || is.list(p_maps))
    stopifnot(identical(p_map_new$labels$title,
                        "Niche Boundary Exceedance"))
    stopifnot(is.null(p_map_new$labels$fill))
    stopifnot(nrow(p_map_new$data) == 8)
    stopifnot(isTRUE(all.equal(p_map_new$data$value, p_map_old$data$value)))
    p_shift <- plot_climniche_map(rf, metric = "niche_distance_change")
    shift_limits <- p_shift$scales$get_scales("fill")$limits
    stopifnot(identical(p_shift$labels$title,
                        "Niche Distance Shift"))
    stopifnot(isTRUE(all.equal(abs(shift_limits[1]), abs(shift_limits[2]))))
    p_custom <- plot_climniche_map(
      rf, metric = "climate_change_amount", limits = c(0, 3),
      extent = c(0, 4, 0, 4), degree_labels = "hemisphere"
    )
    region_boundary <- data.frame(
      x = c(0.5, 3.5, 3.5, 0.5, 0.5),
      y = c(0.5, 0.5, 3.5, 3.5, 0.5),
      group = 1L
    )
    p_region <- plot_climniche_map(
      rf,
      metric = "climate_change_amount",
      study_region = region_boundary
    )
    stopifnot(isTRUE(all.equal(
      p_custom$scales$get_scales("fill")$limits, c(0, 3)
    )))
    stopifnot(inherits(p_region, "ggplot"))
    stopifnot(any(vapply(p_region$layers, function(layer) {
      inherits(layer$geom, "GeomPath")
    }, logical(1))))
    stopifnot(inherits(plot_climniche_variable_contribution(away), "ggplot"))
  }
}

if (requireNamespace("terra", quietly = TRUE)) {
  single_current <- terra::rast(
    nrows = 4, ncols = 4, xmin = 0, xmax = 4,
    ymin = 0, ymax = 4, crs = "EPSG:4326"
  )
  single_values <- seq_len(terra::ncell(single_current))
  single_values[1] <- Inf
  terra::values(single_current) <- single_values
  single_future <- single_current
  terra::values(single_future) <- seq_len(terra::ncell(single_future)) + 0.5
  single_fit <- fit_climniche_terra(single_current, single_future)
  stopifnot(ncol(single_fit$current) == 1L)
  stopifnot(nrow(single_fit$current) == 15L)
  stopifnot(!single_fit$raster_complete[1])

  r1 <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
                    ymin = 0, ymax = 4, crs = "EPSG:4326")
  terra::values(r1) <- seq_len(terra::ncell(r1))
  r2 <- r1
  terra::values(r2) <- rev(seq_len(terra::ncell(r2)))
  cur <- c(r1, r2)
  fut <- c(r1 + 0.5, r2 + 1)
  names(cur) <- c("temperature", "salinity")
  names(fut) <- c("temperature", "salinity")
  occ <- r1
  terra::values(occ) <- c(rep(0, 7), 0.5, rep(0.7, 4), rep(1, 4))
  domain <- r1
  terra::values(domain) <- c(rep(0, 4), rep(1, 12))

  tf <- fit_climniche_terra(cur, fut, occupied = occ,
                            occupied_threshold = 0.5)
  stopifnot(inherits(tf, "climniche_fit"))
  stopifnot(length(tf$occupied) == 8)
  stopifnot(!any(tf$occupied_weight == 0.5))
  stopifnot(isTRUE(all.equal(sort(unique(tf$occupied_weight[tf$occupied])),
                             c(0.7, 1))))
  stopifnot(methods::is(tf$rasters$climate_change_amount, "SpatRaster"))
  stopifnot(!is.null(tf$rasters$climate_reconfiguration))
  stopifnot(!is.null(tf$rasters$niche_boundary_exceedance))
  stopifnot(identical(names(tf$rasters$climate_change_amount),
                      "climate_change_amount"))
  stopifnot(is.null(tf$rasters$classification))
  stopifnot(is.null(tf$class_lookup))
  stopifnot(identical(tf$call[[1]], as.name("fit_climniche_terra")))
  stopifnot(identical(climniche_table(tf, scope = "all")$cell,
                      which(tf$raster_complete)))
  tf_reordered <- fit_climniche_terra(
    cur, fut[[2:1]], occupied = occ, occupied_threshold = 0.5
  )
  stopifnot(isTRUE(all.equal(
    tf$climate_change_amount,
    tf_reordered$climate_change_amount,
    tolerance = 0
  )))
  tf_domain <- fit_climniche_terra(cur, fut, occupied = occ,
                                   occupied_threshold = 0.5,
                                   domain = domain)
  stopifnot(length(tf_domain$climate_change_amount) == 12)
  stopifnot(length(tf_domain$occupied) == 8)
  stopifnot(isTRUE(all.equal(sort(unique(tf_domain$occupied_weight[tf_domain$occupied])),
                             c(0.7, 1))))

  # Future missing values do not alter the current reference centre.
  reference_current <- terra::rast(
    nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2
  )
  terra::values(reference_current) <- c(1, 2, 3, 100)
  reference_future <- reference_current + 1
  reference_future_values <- terra::values(reference_future)[, 1]
  reference_future_values[4] <- NA_real_
  terra::values(reference_future) <- reference_future_values
  reference_weight <- reference_current
  terra::values(reference_weight) <- c(1, 0, 0, 100)
  reference_fit <- fit_climniche_terra(
    reference_current,
    reference_future,
    occupied = reference_weight,
    scale = FALSE,
    preprocess = FALSE
  )
  stopifnot(isTRUE(all.equal(
    unname(reference_fit$center),
    (1 + 100 * 100) / 101,
    tolerance = 1e-12
  )))
  stopifnot(identical(reference_fit$cell_index, 1:3))
  stopifnot(identical(reference_fit$reference_cell_index, 1:4))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_map_new <- plot_climniche_map(tf, metric = "niche_boundary_exceedance",
                                    occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_map_old <- plot_climniche_map(tf, metric = "outside_niche_exceedance",
                                    occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_maps <- plot_climniche_maps(
      tf, occupied = occ, occupied_only = TRUE,
      occupied_threshold = 0.5
    )
    stopifnot(inherits(p_map_new, "ggplot"))
    stopifnot(inherits(p_map_old, "ggplot"))
    stopifnot(inherits(p_maps, "patchwork") || is.list(p_maps))
    stopifnot(nrow(p_map_new$data) == 8)
    stopifnot(isTRUE(all.equal(p_map_new$data$value, p_map_old$data$value)))
  }
}
