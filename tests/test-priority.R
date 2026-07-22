if (!exists("fit_climniche", mode = "function") &&
    requireNamespace("climniche", quietly = TRUE)) {
  library(climniche)
}

priority_internal <- function(name) {
  if (exists(name, mode = "function")) {
    return(get(name, mode = "function"))
  }
  getFromNamespace(name, "climniche")
}

pareto_rank_2d <- priority_internal(".pareto_rank_2d")

# Pareto depth must preserve ties and weak dominance.
first <- c(3, 3, 2, 1, 3, 1)
second <- c(3, 3, 2, 1, 1, 3)
stopifnot(identical(
  pareto_rank_2d(first, second),
  c(1L, 1L, 2L, 3L, 2L, 2L)
))

brute_pareto_rank <- function(first, second) {
  criteria <- cbind(first, second)
  out <- integer(nrow(criteria))
  remaining <- seq_len(nrow(criteria))
  rank <- 1L
  while (length(remaining)) {
    dominated <- vapply(remaining, function(i) {
      candidates <- setdiff(remaining, i)
      any(vapply(candidates, function(j) {
        all(criteria[j, ] >= criteria[i, ]) &&
          any(criteria[j, ] > criteria[i, ])
      }, logical(1)))
    }, logical(1))
    front <- remaining[!dominated]
    out[front] <- rank
    remaining <- setdiff(remaining, front)
    rank <- rank + 1L
  }
  out
}

set.seed(92)
first_test <- sample(1:8, 80, replace = TRUE)
second_test <- sample(1:8, 80, replace = TRUE)
stopifnot(identical(
  pareto_rank_2d(first_test, second_test),
  brute_pareto_rank(first_test, second_test)
))
stopifnot(identical(
  pareto_rank_2d(first_test, second_test),
  pareto_rank_2d(exp(first_test), second_test^3)
))

priority_sim <- simulate_climniche(n = 600, p = 6, seed = 19)
priority_fit <- fit_climniche(
  priority_sim$current,
  priority_sim$future_away,
  occupied = priority_sim$occupied,
  sensitivity = priority_sim$sensitivity,
  preprocess = FALSE
)
priority <- climniche_priority(priority_fit)
stopifnot(inherits(priority, "climniche_priority"))
stopifnot(sum(priority$table$included) == sum(
  priority_fit$occupied_weight > 0 & priority_fit$niche_distance_change > 0
))
stopifnot(all(
  priority$table$niche_distance_change[priority$table$included] >= 0
))
stopifnot(min(priority$table$pareto_rank, na.rm = TRUE) == 1L)
stopifnot(max(priority$table$relative_priority, na.rm = TRUE) == 1)
stopifnot(min(priority$table$relative_priority, na.rm = TRUE) == 0)
rank_order <- order(priority$table$pareto_rank[priority$table$included])
stopifnot(all(diff(
  priority$table$relative_priority[priority$table$included][rank_order]
) <= 0))
stopifnot(priority$front_sizes$n[1L] == sum(
  priority$table$pareto_rank == 1L,
  na.rm = TRUE
))
priority_summary <- summary(priority)
stopifnot(inherits(priority_summary, "summary.climniche_priority"))
stopifnot(isTRUE(all.equal(
  priority_summary$diagnostics$first_front_fraction,
  priority$front_sizes$n[1L] / sum(priority$table$included),
  tolerance = 0
)))
stopifnot(
  priority_summary$diagnostics$unique_exposure_values > 1L,
  priority_summary$diagnostics$unique_criterion_values == 1L
)

all_reference_priority <- climniche_priority(
  priority_fit,
  positive_only = FALSE
)
stopifnot(sum(all_reference_priority$table$included) == sum(
  priority_fit$occupied_weight > 0
))

boundary_priority <- climniche_priority(
  priority_fit,
  exposure = "niche_boundary_exceedance"
)
stopifnot(all(
  boundary_priority$table$niche_boundary_exceedance[
    boundary_priority$table$included
  ] > 0
))

# A user criterion can be minimized without changing the stored values.
management_cost <- seq_len(nrow(priority_fit$current))
cost_priority <- climniche_priority(
  priority_fit,
  exposure = "climate_reconfiguration",
  criterion = management_cost,
  criterion_name = "Management cost",
  criterion_direction = "minimize",
  scope = "all"
)
stopifnot(identical(
  cost_priority$table$decision_criterion,
  as.numeric(management_cost)
))
stopifnot(identical(cost_priority$criterion_direction, "minimize"))
stopifnot(
  summary(cost_priority)$diagnostics$unique_criterion_values > 1L
)
cost_rows <- cost_priority$table$included
stopifnot(identical(
  cost_priority$table$pareto_rank[cost_rows],
  pareto_rank_2d(
    -management_cost[cost_rows],
    priority_fit$climate_reconfiguration[cost_rows]
  )
))

# Legacy exposure names resolve to the primary fields.
legacy_priority <- climniche_priority(
  priority_fit,
  exposure = "composition_change"
)
stopifnot(identical(legacy_priority$exposure, "climate_reconfiguration"))
stopifnot(isTRUE(all.equal(
  legacy_priority$table$climate_reconfiguration,
  priority_fit$climate_reconfiguration,
  tolerance = 0
)))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  priority_plane <- plot_climniche_priority(priority, type = "plane")
  stopifnot(inherits(priority_plane, "ggplot"))
  stopifnot(inherits(plot(priority, type = "plane"), "ggplot"))
}

bad_name <- try(climniche_priority(
  priority_fit,
  criterion_name = character()
), silent = TRUE)
stopifnot(inherits(bad_name, "try-error"))

if (requireNamespace("terra", quietly = TRUE)) {
  set.seed(93)
  template <- terra::rast(
    nrows = 5,
    ncols = 6,
    xmin = -6,
    xmax = 30,
    ymin = 30,
    ymax = 46,
    crs = "EPSG:4326"
  )
  current <- c(template, template, template, template, template, template)
  terra::values(current) <- matrix(rnorm(terra::ncell(template) * 6), ncol = 6)
  names(current) <- paste0("climate_", seq_len(6))
  future <- current
  terra::values(future) <- terra::values(current) + matrix(
    rep(c(0.5, 0.2, -0.1, 0.3, -0.2, 0.4), terra::ncell(template)),
    ncol = 6,
    byrow = TRUE
  )
  suitability <- template
  suitability_values <- seq(0.05, 0.95, length.out = terra::ncell(template))
  terra::values(suitability) <- suitability_values

  spatial_fit <- fit_climniche_terra(
    current,
    future,
    occupied = suitability,
    occupied_threshold = 0.3,
    sensitivity = rep(1, 6),
    preprocess = FALSE
  )
  spatial_priority <- climniche_priority(spatial_fit)
  stopifnot(inherits(spatial_priority$rasters$pareto_rank, "SpatRaster"))
  rank_values <- terra::values(spatial_priority$rasters$pareto_rank)[, 1L]
  expected_ranked <- suitability_values > 0.3 &
    spatial_fit$niche_distance_change > 0
  stopifnot(all(is.finite(rank_values[expected_ranked])))
  stopifnot(all(is.na(rank_values[!expected_ranked])))

  spatial_value_priority <- climniche_priority(
    spatial_fit,
    criterion = suitability,
    criterion_name = "Habitat value",
    scope = "all"
  )
  stopifnot(isTRUE(all.equal(
    spatial_value_priority$table$decision_criterion,
    suitability_values,
    tolerance = 0
  )))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    priority_map <- plot_climniche_priority(
      spatial_priority,
      type = "map",
      degree_labels = "hemisphere"
    )
    stopifnot(inherits(priority_map, "ggplot"))
  }
}
