if (!exists("fit_climniche", mode = "function") &&
    requireNamespace("climniche", quietly = TRUE)) {
  library(climniche)
}

normalise_class_test <- if (exists(".normalise_class", mode = "function")) {
  .normalise_class
} else {
  get(".normalise_class", envir = asNamespace("climniche"))
}

sim <- simulate_climniche(n = 500, seed = 10)

toward <- fit_climniche(
  current = sim$current,
  future = sim$future_toward,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity
)

away <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity
)

away_numeric_occupied <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = as.numeric(sim$occupied),
  occupied_threshold = 0.5,
  center = sim$center,
  sensitivity = sim$sensitivity
)

stopifnot(inherits(toward, "climniche_fit"))
stopifnot(identical(away$occupied, away_numeric_occupied$occupied))
stopifnot(length(toward$climate_change_amount) == nrow(sim$current))
stopifnot(length(toward$climate_reconfiguration) == nrow(sim$current))
stopifnot(length(toward$niche_boundary_exceedance) == nrow(sim$current))
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
stopifnot(all(c("occupied_weight", "climate_reconfiguration",
                "niche_boundary_exceedance", "composition_change",
                "outside_niche_exceedance", "radial_direction",
                "boundary_status") %in%
                names(climniche_table(away))))
stopifnot(identical(levels(away$radial_direction), c(
  "Toward realised niche centre",
  "Limited Niche Distance Shift",
  "Away from realised niche centre"
)))
stopifnot(identical(levels(away$boundary_status), c(
  "Within empirical niche boundary",
  "Beyond empirical niche boundary"
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
stopifnot(isTRUE(all.equal(weighted_fit$center, expected_center)))
stopifnot(!isTRUE(all.equal(weighted_fit$center, colMeans(simple_current[1:2, ]))))

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
stopifnot(isTRUE(all.equal(length_n_weight_fit$center,
                           expected_length_n_center)))

negative_weight <- try(fit_climniche(
  simple_current,
  simple_future,
  occupied = c(1, -1, 0, 0),
  sensitivity = c(1, 1),
  scale = FALSE
), silent = TRUE)
stopifnot(inherits(negative_weight, "try-error"))

custom_class <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity,
  tolerance = 1e6,
  stable_climate_change = 1e6,
  stable_reconfiguration = 1e6,
  boundary_exceedance_tolerance = 1e6
)
stopifnot(!identical(as.character(away$classification),
                     as.character(custom_class$classification)))
stopifnot(all(custom_class$classification == "Limited niche relative change"))
stopifnot(all(custom_class$radial_direction == "Limited Niche Distance Shift"))
stopifnot(all(custom_class$boundary_status ==
                "Within empirical niche boundary"))
stopifnot(as.character(normalise_class_test("Limited climate niche change")) ==
            "Limited niche relative change")
stopifnot(as.character(normalise_class_test("little climate niche change")) ==
            "Limited niche relative change")
stopifnot(isTRUE(all.equal(custom_class$classification_settings$tolerance,
                           1e6)))
stopifnot(isTRUE(all.equal(
  custom_class$classification_settings$stable_climate_change,
  1e6
)))
stopifnot(isTRUE(all.equal(
  custom_class$classification_settings$stable_reconfiguration,
  1e6
)))
stopifnot(isTRUE(all.equal(
  custom_class$classification_settings$boundary_exceedance_tolerance,
  1e6
)))
stopifnot(is.finite(away$classification_settings$tolerance))
stopifnot(is.finite(away$classification_settings$stable_climate_change))
stopifnot(is.finite(away$classification_settings$stable_reconfiguration))
limited_summary <- climniche_summary(custom_class)
stopifnot(isTRUE(all.equal(limited_summary$prop_stable, 1)))
stopifnot(isTRUE(all.equal(limited_summary$prop_niche_convergence, 0)))
stopifnot(isTRUE(all.equal(limited_summary$prop_niche_divergence, 0)))
stopifnot(isTRUE(all.equal(limited_summary$prop_niche_exceedance, 0)))

A_full <- matrix(c(1, 0.25, 0.25, 1.4), 2, 2)
fit_full <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  A = A_full
)
stopifnot(isTRUE(all.equal(rowSums(fit_full$variable_contribution),
                           fit_full$psi_future - fit_full$psi_current,
                           tolerance = 1e-7)))
bad_A <- matrix(c(1, 2, 0, 1), 2, 2)
bad_fit <- try(fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  A = bad_A
), silent = TRUE)
stopifnot(inherits(bad_fit, "try-error"))
stopifnot(is.logical(away$mixed_variable_response))
stopifnot(!("mixed variable response" %in% levels(away$classification)))
stopifnot("Climatic Reconfiguration with limited Niche Distance Shift" %in%
            levels(away$classification))

diagram <- climniche_diagram_data(away, max_arrows = 50)
stopifnot(inherits(diagram, "climniche_diagram_data"))
stopifnot(nrow(diagram$cells) == length(away$occupied))
stopifnot(nrow(diagram$arrows) <= 50)
stopifnot(all(c("current_axis1", "future_axis1", "class",
                "climate_reconfiguration", "niche_boundary_exceedance",
                "occupied_weight") %in% names(diagram$cells)))
showcase <- climniche_showcase_data(away, max_points = 80)
summary_data <- climniche_summary_figure_data(away, max_points = 80)
stopifnot(inherits(showcase, "climniche_showcase_data"))
stopifnot(inherits(summary_data, "climniche_showcase_data"))
stopifnot(nrow(showcase$plane) <= 80)
stopifnot(all(c("x_mid", "y_mid", "cell_weight", "total_weight") %in%
                names(showcase$plane_bins)))
stopifnot(all(c("x_mid", "y_mid", "total_weight") %in%
                names(showcase$reconfiguration_bins)))
stopifnot(all(c("metric", "value", "weight") %in% names(showcase$metrics)))
stopifnot(all(c("proportion", "count", "weight") %in%
                names(showcase$classes)))
stopifnot(nrow(showcase$classes) == 5L)
stopifnot(all(c("descriptor", "level", "proportion") %in%
                names(showcase$descriptors)))
stopifnot(nrow(showcase$descriptors) == 5L)
stopifnot(all(c("metric", "x", "width", "proportion") %in%
                names(showcase$metric_histograms)))
stopifnot("mean_absolute_share" %in% names(showcase$variables))
stopifnot(all(c("boundary_quantile", "prop_exceeded") %in%
                names(showcase$boundary)))
report_all_classes <- climniche_report(away, scope = "current")
formal_metric_names <- c(
  "Climatic Displacement",
  "Niche Distance Shift",
  "Climatic Reconfiguration",
  "Niche Boundary Exceedance"
)
stopifnot(all(formal_metric_names %in%
                report_all_classes$metric_definitions$metric))
stopifnot(all(formal_metric_names %in%
                report_all_classes$metric_summary$metric))
tmp_report <- tempfile(fileext = ".md")
write_climniche_report(report_all_classes, tmp_report)
report_text <- paste(readLines(tmp_report, warn = FALSE), collapse = "\n")
stopifnot(all(vapply(formal_metric_names, grepl, logical(1),
                     x = report_text, fixed = TRUE)))
stopifnot(nrow(report_all_classes$descriptor_summary) == 5L)
stopifnot(nrow(report_all_classes$class_summary) == 5L)
stopifnot(all(as.character(report_all_classes$class_summary$class) ==
                c(
                  "Limited niche relative change",
                  "Closer to current niche",
                  "Farther from current niche",
                  "Outside current niche boundary",
                  paste(
                    "Climatic Reconfiguration with limited",
                    "Niche Distance Shift"
                  )
                )))

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
  p_exposure_classes <- plot_climniche_exposure(
    away, colour_by = "classification"
  )
  stopifnot(inherits(p_exposure, "ggplot"))
  stopifnot(inherits(p_exposure_classes, "ggplot"))
  stopifnot(identical(p_exposure$labels$x, "Climatic Displacement"))
  stopifnot(identical(p_exposure$labels$y, "Niche Distance Shift"))
  showcase_plot <- plot_climniche_showcase(away, max_points = 80)
  summary_figure_plot <- plot_climniche_summary_figure(away, max_points = 80)
  stopifnot(inherits(showcase_plot, "patchwork") || is.list(showcase_plot))
  stopifnot(inherits(summary_figure_plot, "patchwork") ||
              is.list(summary_figure_plot))
}

if (requireNamespace("raster", quietly = TRUE)) {
  r1 <- raster::raster(nrows = 4, ncols = 4, xmn = 0, xmx = 4,
                       ymn = 0, ymx = 4, crs = "+proj=longlat +datum=WGS84")
  raster::values(r1) <- seq_len(raster::ncell(r1))
  r2 <- r1
  raster::values(r2) <- rev(seq_len(raster::ncell(r2)))
  cur <- raster::stack(r1, r2)
  fut <- raster::stack(r1 + 0.5, r2 + 1)
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
  rf_domain <- fit_climniche_raster(cur, fut, occupied = occ,
                                    occupied_threshold = 0.5,
                                    domain = domain)
  stopifnot(length(rf_domain$climate_change_amount) == 12)
  stopifnot(length(rf_domain$occupied) == 8)
  stopifnot(isTRUE(all.equal(sort(unique(rf_domain$occupied_weight[rf_domain$occupied])),
                             c(0.7, 1))))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_map_new <- plot_climniche_map(rf, metric = "niche_boundary_exceedance",
                                    occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_map_old <- plot_climniche_map(rf, metric = "outside_niche_exceedance",
                                    occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_cls <- plot_climniche_classes(rf, occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_cls_all <- plot_climniche_classes(
      rf, occupied = occ, occupied_only = TRUE,
      occupied_threshold = 0.5, class_display = "all"
    )
    p_maps <- plot_climniche_maps(
      rf, occupied = occ, occupied_only = TRUE,
      occupied_threshold = 0.5, degree_labels = "hemisphere"
    )
    stopifnot(inherits(p_map_new, "ggplot"))
    stopifnot(inherits(p_map_old, "ggplot"))
    stopifnot(inherits(p_cls, "ggplot"))
    stopifnot(inherits(p_cls_all, "ggplot"))
    stopifnot(inherits(p_maps, "patchwork") || is.list(p_maps))
    stopifnot(identical(p_map_new$labels$title,
                        "Niche Boundary Exceedance"))
    stopifnot(nrow(p_map_new$data) == 8)
    stopifnot(nrow(p_cls$data) == 8)
    stopifnot(isTRUE(all.equal(p_map_new$data$value, p_map_old$data$value)))
    stopifnot(isTRUE(all.equal(
      p_map_new$scales$get_scales("fill")$limits[1], 0
    )))
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
    if (requireNamespace("sf", quietly = TRUE)) {
      region_sf <- sf::st_sf(
        geometry = sf::st_sfc(
          sf::st_polygon(list(as.matrix(region_boundary[, c("x", "y")]))),
          crs = 4326
        )
      )
      p_region_sf <- plot_climniche_map(
        rf,
        metric = "climate_change_amount",
        study_region = region_sf
      )
      stopifnot(inherits(p_region_sf, "ggplot"))
    }
    stopifnot(inherits(plot_climniche_variable_contribution(away), "ggplot"))
  }
}

if (requireNamespace("terra", quietly = TRUE)) {
  r1 <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
                    ymin = 0, ymax = 4, crs = "EPSG:4326")
  terra::values(r1) <- seq_len(terra::ncell(r1))
  r2 <- r1
  terra::values(r2) <- rev(seq_len(terra::ncell(r2)))
  cur <- c(r1, r2)
  fut <- c(r1 + 0.5, r2 + 1)
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
  tf_domain <- fit_climniche_terra(cur, fut, occupied = occ,
                                   occupied_threshold = 0.5,
                                   domain = domain)
  stopifnot(length(tf_domain$climate_change_amount) == 12)
  stopifnot(length(tf_domain$occupied) == 8)
  stopifnot(isTRUE(all.equal(sort(unique(tf_domain$occupied_weight[tf_domain$occupied])),
                             c(0.7, 1))))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_map_new <- plot_climniche_map(tf, metric = "niche_boundary_exceedance",
                                    occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_map_old <- plot_climniche_map(tf, metric = "outside_niche_exceedance",
                                    occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_cls <- plot_climniche_classes(tf, occupied = occ, occupied_only = TRUE,
                                    occupied_threshold = 0.5)
    p_maps <- plot_climniche_maps(
      tf, occupied = occ, occupied_only = TRUE,
      occupied_threshold = 0.5
    )
    stopifnot(inherits(p_map_new, "ggplot"))
    stopifnot(inherits(p_map_old, "ggplot"))
    stopifnot(inherits(p_cls, "ggplot"))
    stopifnot(inherits(p_maps, "patchwork") || is.list(p_maps))
    stopifnot(nrow(p_map_new$data) == 8)
    stopifnot(nrow(p_cls$data) == 8)
    stopifnot(isTRUE(all.equal(p_map_new$data$value, p_map_old$data$value)))
  }
}
