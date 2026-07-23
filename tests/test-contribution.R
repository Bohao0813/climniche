if (!exists("fit_climniche", mode = "function") &&
    requireNamespace("climniche", quietly = TRUE)) {
  library(climniche)
}

current <- matrix(0, nrow = 5, ncol = 2)
future <- rbind(
  c(2, 1),
  c(1, 3),
  c(-4, 1),
  c(1, 1),
  c(0, 0)
)
colnames(current) <- colnames(future) <- c("temperature", "water_balance")
reference_weight <- c(1, 0.5, 0, 1, 1)

fit <- fit_climniche(
  current,
  future,
  occupied = reference_weight,
  center = c(temperature = 0, water_balance = 0),
  A = diag(2),
  scale = FALSE,
  preprocess = FALSE,
  tolerance = 0
)
contribution <- climniche_dominant_contribution(fit)
stopifnot(inherits(contribution, "climniche_contribution"))
stopifnot(identical(
  contribution$table$dominant_variable,
  c("temperature", "water_balance", NA, "Tied", NA)
))
stopifnot(isTRUE(all.equal(
  contribution$table$dominant_share[1:2],
  c(4 / 5, 9 / 10),
  tolerance = 1e-12
)))
stopifnot(isTRUE(contribution$table$tied[4]))
stopifnot(is.na(contribution$table$dominant_contribution[4]))
stopifnot(!contribution$table$included[3])
stopifnot(isTRUE(all.equal(
  rowSums(fit$variable_contribution),
  fit$psi_future - fit$psi_current,
  tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  rowSums(fit$variable_contribution),
  fit$niche_distance_change *
    (fit$niche_radius_future + fit$niche_radius_current),
  tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  sum(contribution$summary$mean_absolute_share),
  1,
  tolerance = 1e-12
)))
expected_signed <- c(5.5, 6.5) / 3.5
stopifnot(isTRUE(all.equal(
  contribution$summary$mean_signed_contribution[
    match(c("temperature", "water_balance"), contribution$summary$variable)
  ],
  expected_signed,
  tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  contribution$summary$positive_contribution_fraction[
    match(c("temperature", "water_balance"), contribution$summary$variable)
  ],
  rep(2.5 / 3.5, 2),
  tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  contribution$summary$dominant_weight_fraction[
    match(c("temperature", "water_balance"), contribution$summary$variable)
  ],
  c(1 / 2.5, 0.5 / 2.5),
  tolerance = 1e-12
)))
stopifnot(identical(summary(contribution), contribution$summary))
figure_data <- climniche_summary_figure_data(
  fit,
  scope = "current",
  top_variables = 2
)
summary_share <- contribution$summary$mean_absolute_share[
  match(
    as.character(figure_data$variables$variable),
    contribution$summary$variable
  )
]
stopifnot(isTRUE(all.equal(
  summary_share,
  figure_data$variables$mean_absolute_share,
  tolerance = 1e-12
)))
report <- climniche_report(fit, top_variables = 2)
report_share <- report$top_variables$mean_absolute_share[
  match(
    contribution$summary$variable,
    report$top_variables$variable
  )
]
stopifnot(isTRUE(all.equal(
  report_share,
  contribution$summary$mean_absolute_share,
  tolerance = 1e-12
)))
stopifnot(all(c(
  "mean_contribution", "abs_mean_contribution"
) %in% names(report$top_variables)))
stopifnot(inherits(print(contribution), "climniche_contribution"))

# Cross-variable terms remain part of the exact contribution identity.
A <- matrix(c(1, 0.35, 0.35, 1.2), nrow = 2)
fit_full <- fit_climniche(
  current,
  future,
  occupied = rep(1, 5),
  center = c(temperature = 0, water_balance = 0),
  A = A,
  scale = FALSE,
  preprocess = FALSE,
  tolerance = 0
)
full_contribution <- climniche_dominant_contribution(fit_full, scope = "all")
stopifnot(isTRUE(all.equal(
  full_contribution$table$squared_niche_distance_change,
  fit_full$psi_future - fit_full$psi_current,
  tolerance = 1e-12
)))
stopifnot(identical(
  full_contribution$table$squared_niche_distance_change,
  full_contribution$table$niche_potential_change
))

zero_fit <- fit_climniche(
  current,
  current,
  occupied = rep(1, 5),
  center = c(temperature = 0, water_balance = 0),
  A = diag(2),
  scale = FALSE,
  preprocess = FALSE,
  tolerance = 0
)
zero_contribution <- climniche_dominant_contribution(zero_fit)
stopifnot(all(is.na(zero_contribution$table$dominant_variable)))
stopifnot(all(is.na(zero_contribution$summary$mean_absolute_share)))

if (requireNamespace("terra", quietly = TRUE)) {
  template <- terra::rast(
    nrows = 2,
    ncols = 3,
    xmin = -5,
    xmax = 25,
    ymin = 30,
    ymax = 45,
    crs = "EPSG:4326"
  )
  current_raster <- c(template, template)
  terra::values(current_raster) <- matrix(0, nrow = 6, ncol = 2)
  names(current_raster) <- c("temperature", "water_balance")
  future_raster <- current_raster
  terra::values(future_raster) <- rbind(
    c(2, 1), c(1, 3), c(4, 1),
    c(1, 2), c(3, 1), c(1, 4)
  )
  reference <- template
  terra::values(reference) <- c(1, 0.7, 0, 0.4, 0.8, 0)

  spatial_fit <- fit_climniche_terra(
    current_raster,
    future_raster,
    occupied = reference,
    center = c(temperature = 0, water_balance = 0),
    A = diag(2),
    scale = FALSE,
    preprocess = FALSE,
    tolerance = 0
  )
  spatial_contribution <- climniche_dominant_contribution(spatial_fit)
  stopifnot(inherits(
    spatial_contribution$rasters$dominant_variable,
    "SpatRaster"
  ))
  dominant_values <- terra::values(
    spatial_contribution$rasters$dominant_variable
  )[, 1L]
  stopifnot(all(is.na(dominant_values[c(3, 6)])))
  stopifnot(all(is.finite(dominant_values[c(1, 2, 4, 5)])))
  stopifnot(identical(
    terra::values(
      spatial_contribution$rasters$squared_niche_distance_change
    )[, 1L],
    terra::values(
      spatial_contribution$rasters$niche_potential_change
    )[, 1L]
  ))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    variable_map <- plot_climniche_dominant_contribution(
      spatial_contribution,
      type = "variable",
      degree_labels = "hemisphere"
    )
    share_map <- plot_climniche_dominant_contribution(
      spatial_contribution,
      type = "share",
      degree_labels = "hemisphere"
    )
    stopifnot(inherits(variable_map, "ggplot"))
    stopifnot(inherits(share_map, "ggplot"))
    stopifnot(inherits(plot(spatial_contribution, type = "share"), "ggplot"))
    custom_map <- plot_climniche_dominant_contribution(
      spatial_contribution,
      type = "variable",
      colours = c("#0072B2", "#D55E00")
    )
    stopifnot(inherits(custom_map, "ggplot"))
    complete_legend_map <- plot_climniche_dominant_contribution(
      spatial_contribution,
      type = "variable",
      variable_labels = c(
        temperature = "Temperature",
        water_balance = "Water balance"
      ),
      legend_variables = c("temperature", "water_balance")
    )
    fill_scale <- complete_legend_map$scales$get_scales("fill")
    stopifnot(identical(fill_scale$breaks, c("Temperature", "Water balance")))
    stopifnot(identical(fill_scale$drop, FALSE))
  }
}
