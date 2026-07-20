if (!exists("fit_climniche", mode = "function") &&
    requireNamespace("climniche", quietly = TRUE)) {
  library(climniche)
}

current <- matrix(c(
  0, 1, 4,
  1, 3, 2,
  2, 2, 5,
  3, 5, 1,
  4, 4, 3,
  5, 7, 6
), ncol = 3, byrow = TRUE)
colnames(current) <- c("temperature", "salinity", "oxygen")
rownames(current) <- paste0("cell_", seq_len(nrow(current)))
future <- current + matrix(c(0.3, -0.2, 0.1), nrow(current), 3,
                           byrow = TRUE)
weights <- c(0.2, 0.8, 1.0, 0.5, 0, 0.3)
sensitivity <- c(temperature = 2, salinity = 1, oxygen = 0.5)

reference_fit <- fit_climniche(
  current,
  future,
  occupied = weights,
  sensitivity = sensitivity,
  scale = FALSE,
  preprocess = FALSE
)
stopifnot(identical(reference_fit$metric_type, "diag"))
stopifnot(identical(names(reference_fit$center), colnames(current)))
stopifnot(identical(rownames(reference_fit$A), colnames(current)))
stopifnot(identical(colnames(reference_fit$A), colnames(current)))
stopifnot(identical(names(reference_fit$sensitivity_weights),
                    colnames(current)))
stopifnot(identical(names(reference_fit$standardization$mean),
                    colnames(current)))
stopifnot(identical(reference_fit$call[[1]], as.name("fit_climniche")))

# Named rows, columns and sensitivity weights are matched before fitting.
reordered_fit <- fit_climniche(
  current,
  future[rev(rownames(future)), rev(colnames(future)), drop = FALSE],
  occupied = weights,
  sensitivity = sensitivity[rev(names(sensitivity))],
  scale = FALSE,
  preprocess = FALSE
)
stopifnot(isTRUE(all.equal(
  reference_fit$climate_change_amount,
  reordered_fit$climate_change_amount,
  tolerance = 0
)))
stopifnot(isTRUE(all.equal(
  reference_fit$niche_distance_change,
  reordered_fit$niche_distance_change,
  tolerance = 0
)))

bad_variables <- future
colnames(bad_variables)[3] <- "ph"
bad_variables_fit <- try(fit_climniche(
  current,
  bad_variables,
  occupied = weights,
  scale = FALSE,
  preprocess = FALSE
), silent = TRUE)
stopifnot(inherits(bad_variables_fit, "try-error"))

bad_rows <- future
rownames(bad_rows)[6] <- "different_cell"
bad_rows_fit <- try(fit_climniche(
  current,
  bad_rows,
  occupied = weights,
  scale = FALSE,
  preprocess = FALSE
), silent = TRUE)
stopifnot(inherits(bad_rows_fit, "try-error"))

# Named centres and metric matrices follow climate-variable order.
center <- c(temperature = 1.2, salinity = 2.8, oxygen = 3.5)
A <- matrix(c(
  1.5, 0.2, 0.0,
  0.2, 1.0, 0.1,
  0.0, 0.1, 0.8
), nrow = 3, byrow = TRUE,
  dimnames = list(names(center), names(center)))
named_fit <- fit_climniche(
  current,
  future,
  occupied = weights,
  center = center,
  A = A,
  scale = FALSE,
  preprocess = FALSE
)
reverse_names <- rev(names(center))
reordered_named_fit <- fit_climniche(
  current,
  future,
  occupied = weights,
  center = center[reverse_names],
  A = A[reverse_names, reverse_names],
  scale = FALSE,
  preprocess = FALSE
)
stopifnot(isTRUE(all.equal(
  named_fit$niche_boundary_exceedance,
  reordered_named_fit$niche_boundary_exceedance,
  tolerance = 0
)))
direct_potential <- niche_potential(
  current,
  center = center[reverse_names],
  A = A[reverse_names, reverse_names]
)
stopifnot(isTRUE(all.equal(
  direct_potential,
  niche_potential(current, center = center, A = A),
  tolerance = 0
)))

# Named reference weights follow current row names.
named_weights <- stats::setNames(rev(weights), rev(rownames(current)))
named_weight_fit <- fit_climniche(
  current,
  future,
  occupied = named_weights,
  sensitivity = sensitivity,
  scale = FALSE,
  preprocess = FALSE
)
stopifnot(isTRUE(all.equal(
  reference_fit$center,
  named_weight_fit$center,
  tolerance = 0
)))

bad_weight_names <- stats::setNames(weights, paste0("row_", seq_along(weights)))
bad_named_weight_fit <- try(fit_climniche(
  current,
  future,
  occupied = bad_weight_names,
  sensitivity = sensitivity,
  scale = FALSE,
  preprocess = FALSE
), silent = TRUE)
stopifnot(inherits(bad_named_weight_fit, "try-error"))

non_numeric_center <- try(fit_climniche(
  current,
  future,
  occupied = weights,
  center = factor(c("1", "2", "3")),
  A = A,
  scale = FALSE,
  preprocess = FALSE
), silent = TRUE)
stopifnot(inherits(non_numeric_center, "try-error"))

duplicate_current <- current
colnames(duplicate_current) <- c("temperature", "temperature", "oxygen")
unnamed_future <- future
colnames(unnamed_future) <- NULL
duplicate_variables <- try(fit_climniche(
  duplicate_current,
  unnamed_future,
  occupied = weights,
  scale = FALSE,
  preprocess = FALSE
), silent = TRUE)
stopifnot(inherits(duplicate_variables, "try-error"))

# Correlation filtering is invariant to named column order.
correlated_current <- cbind(
  zeta = seq_len(8),
  alpha = 2 * seq_len(8),
  beta = rep(c(1, 2), 4)
)
correlated_future <- correlated_current + 0.1
correlated_fit <- fit_climniche(
  correlated_current,
  correlated_future,
  occupied = rep(1, 8)
)
column_order <- c("beta", "alpha", "zeta")
reordered_correlated_fit <- fit_climniche(
  correlated_current[, column_order],
  correlated_future[, column_order],
  occupied = rep(1, 8)
)
stopifnot(identical(
  sort(correlated_fit$preprocessing$retained_variables),
  sort(reordered_correlated_fit$preprocessing$retained_variables)
))
stopifnot(isTRUE(all.equal(
  correlated_fit$climate_change_amount,
  reordered_correlated_fit$climate_change_amount,
  tolerance = 1e-12
)))

# Binary and continuous weights use the same inverse weighted empirical quantile.
quantile_current <- matrix(seq_len(10), ncol = 1,
                           dimnames = list(NULL, "temperature"))
quantile_fit <- fit_climniche(
  quantile_current,
  quantile_current,
  occupied = rep(2, 10),
  center = 0,
  sensitivity = 1,
  boundary = 0.5,
  scale = FALSE,
  preprocess = FALSE
)
stopifnot(isTRUE(all.equal(quantile_fit$boundary_distance, 5)))
stopifnot(isTRUE(all.equal(
  boundary_exceedance(
    quantile_fit$psi_future,
    quantile_fit$boundary_potential
  ),
  quantile_fit$niche_boundary_exceedance,
  tolerance = 0
)))

factor_coordinates <- diag(3)
dimnames(factor_coordinates) <- list(
  names(sensitivity),
  c("Marg", "Spec1", "Spec2")
)
factor_metric <- niche_metric(
  cnfa = list(
    co = factor_coordinates,
    eig = c(Marg = 2, Spec1 = 1, Spec2 = 0.5)
  ),
  type = "factor"
)
stopifnot(identical(dim(factor_metric), c(3L, 3L)))
stopifnot(identical(rownames(factor_metric), names(sensitivity)))

# A requested factor metric must not change method when CNFA data are missing.
missing_factor <- try(niche_metric(
  sensitivity = sensitivity,
  type = "factor"
), silent = TRUE)
stopifnot(inherits(missing_factor, "try-error"))

bad_scale <- try(fit_climniche(
  current,
  future,
  occupied = weights,
  scale = NA,
  preprocess = FALSE
), silent = TRUE)
stopifnot(inherits(bad_scale, "try-error"))

bad_sd <- try(fit_climniche(
  current,
  future,
  occupied = weights,
  global_sd = c(1, NA, 1),
  preprocess = FALSE
), silent = TRUE)
stopifnot(inherits(bad_sd, "try-error"))

bad_report_count <- try(climniche_report(reference_fit, top_variables = 0),
                        silent = TRUE)
stopifnot(inherits(bad_report_count, "try-error"))

bad_figure_quantiles <- try(climniche_summary_figure_data(
  reference_fit,
  boundary_probs = c(0, 0.95)
), silent = TRUE)
stopifnot(inherits(bad_figure_quantiles, "try-error"))

bad_seed <- try(simulate_climniche(seed = 1.5), silent = TRUE)
stopifnot(inherits(bad_seed, "try-error"))

current_report <- climniche_report(reference_fit, scope = "current")
all_report <- climniche_report(reference_fit, scope = "all")
stopifnot(grepl("total reference weight", current_report$interpretation[2],
                fixed = TRUE))
stopifnot(grepl("analysed cells", all_report$interpretation[2], fixed = TRUE))
stopifnot(identical(current_report$settings$metric_type, "diag"))
stopifnot(isTRUE(all.equal(
  current_report$metric_weights$sensitivity_weight,
  reference_fit$sensitivity_weights,
  check.attributes = FALSE
)))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  weighted_plot <- plot_climniche_distribution(
    reference_fit,
    metric = "climate_change_amount",
    scope = "current"
  )
  histogram <- ggplot2::ggplot_build(weighted_plot)$data[[1]]
  stopifnot(isTRUE(all.equal(
    sum(histogram$count),
    sum(reference_fit$occupied_weight[reference_fit$occupied])
  )))
  stopifnot(identical(weighted_plot$labels$y, "Weighted count"))
}
