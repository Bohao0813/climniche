#' Screen ecological criteria against climate exposure
#'
#' @param x A fitted `climniche_fit` object.
#' @param exposure Climatic quantity used as the exposure objective. Available
#'   choices are `"niche_distance_change"`, `"niche_boundary_exceedance"`,
#'   `"climate_reconfiguration"`, and `"climate_change_amount"`. Legacy metric
#'   aliases are accepted.
#' @param criterion Optional second decision criterion. A numeric vector may
#'   contain one value per evaluated row or, for a spatial fit, one value per
#'   raster cell. A matching one-layer RasterLayer or SpatRaster is also
#'   accepted. When omitted, current reference weights are used.
#' @param criterion_name Display name for `criterion`.
#' @param criterion_direction Whether larger or smaller criterion values are
#'   preferred.
#' @param scope `"current"` ranks cells with positive current reference weight;
#'   `"all"` ranks every evaluated cell with finite criteria.
#' @param positive_only If `TRUE`, only cells with a positive value of the
#'   selected exposure quantity are ranked. Set to `FALSE` to retain zero and
#'   negative values.
#' @param exposure_direction Whether larger or smaller exposure values are
#'   preferred in the Pareto comparison. Use `"maximize"` to screen for
#'   greater exposure and `"minimize"` to screen for lower exposure.
#'
#' @return A `climniche_priority` object containing the two decision criteria,
#'   Pareto ranks, Pareto depth score (`pareto_depth_score`) and, for spatial
#'   fits, map layers. `relative_priority` is retained as a compatibility alias
#'   for `pareto_depth_score`.
#'
#' @details
#' The function ranks cells on two objectives: the selected climatic exposure
#' quantity and one reference or decision criterion. Cell \eqn{i} dominates cell
#' \eqn{k} when it is at least as preferred on both objectives and strictly
#' preferred on one. Non-dominated cells form Pareto rank 1. Removing that front
#' and repeating the comparison produces ranks 2, 3, and subsequent fronts.
#'
#' If an analysis contains \eqn{K} Pareto fronts, `pareto_depth_score` is
#' \eqn{(K-r_i)/(K-1)} for the rank \eqn{r_i} of cell \eqn{i}; it is 1 when all
#' cells occupy one front. The score does not order cells within a front and is
#' not comparable among separate analyses. Pareto dominance is invariant to
#' monotonic rescaling, so the two criteria are not combined with fitted
#' weights.
#'
#' When current reference weights are constant, the default second criterion
#' does not distinguish cells and ranking is determined by exposure alone.
#'
#' Only one exposure quantity is used at a time. This avoids counting Climatic
#' Displacement, Niche Distance Shift and Climatic Reconfiguration as independent
#' objectives even though their values satisfy the fitted geometric identity.
#' When the second criterion represents ecological value, maximising a
#' positive Niche Distance Shift identifies cells where high ecological value
#' coincides with movement away from the realised niche centre. Minimising
#' Climatic Displacement identifies cells where high ecological value
#' coincides with less local climatic change. The two directions answer
#' different screening questions; neither is a complete conservation ranking.
#'
#' `summary()` reports the first-front fraction and the Spearman correlation
#' between the preference-oriented forms of the two objectives. These describe
#' how strongly the objectives separate the ranked cells; they are not
#' inferential tests.
#'
#' @references Tracey JA, Rochester CJ, Hathaway SA, et al. (2018).
#'   Prioritizing conserved areas threatened by wildfire and fragmentation for
#'   monitoring and management. *PLOS ONE*, 13, e0200203.
#'   \doi{10.1371/journal.pone.0200203}
#'
#'   Sacre E, Bode M, Weeks R, Pressey RL (2019). The context dependence of
#'   frontier versus wilderness conservation priorities. *Conservation
#'   Letters*, 12, e12632. \doi{10.1111/conl.12632}
#'
#' @examples
#' sim <- simulate_climniche(n = 500, p = 6, seed = 18)
#' fit <- fit_climniche(
#'   sim[["current"]],
#'   sim[["future_away"]],
#'   occupied = sim[["occupied"]],
#'   sensitivity = sim[["sensitivity"]]
#' )
#' priority <- climniche_priority(fit)
#' priority
#' priority_table <- priority[["table"]]
#' head(priority_table[priority_table[["included"]], ])
#' summary(priority)
#' @export
climniche_priority <- function(
    x,
    exposure = c(
      "niche_distance_change",
      "niche_boundary_exceedance",
      "climate_reconfiguration",
      "climate_change_amount",
      "outside_niche_exceedance",
      "composition_change"
    ),
    criterion = NULL,
    criterion_name = NULL,
    criterion_direction = c("maximize", "minimize"),
    scope = c("current", "all"),
    positive_only = TRUE,
    exposure_direction = c("maximize", "minimize")) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  exposure <- match.arg(exposure)
  criterion_direction <- match.arg(criterion_direction)
  exposure_direction <- match.arg(exposure_direction)
  scope <- match.arg(scope)
  positive_only <- .check_flag(positive_only, "positive_only")
  if (!is.null(criterion_name) &&
      (!is.character(criterion_name) || length(criterion_name) != 1L ||
       is.na(criterion_name) || !nzchar(criterion_name))) {
    stop("criterion_name must be one non-empty character value.",
         call. = FALSE)
  }

  exposure_data <- .priority_exposure_values(
    x,
    exposure = exposure,
    positive_only = positive_only
  )
  criterion_data <- .priority_criterion_values(
    x,
    criterion = criterion,
    criterion_name = criterion_name
  )
  reference_weight <- .fit_reference_weights(x)
  n <- length(exposure_data$values)
  if (length(reference_weight) != n) {
    stop("x contains incompatible reference weights.", call. = FALSE)
  }

  included <- is.finite(exposure_data$values) &
    is.finite(criterion_data$values)
  if (positive_only) {
    included <- included & exposure_data$values > 0
  }
  if (identical(scope, "current")) {
    included <- included & reference_weight > 0
  }
  if (!any(included)) {
    stop("No cells have finite screening criteria in the selected scope.",
         call. = FALSE)
  }

  criterion_preference <- criterion_data$values[included]
  if (identical(criterion_direction, "minimize")) {
    criterion_preference <- -criterion_preference
  }
  exposure_preference <- exposure_data$values[included]
  if (identical(exposure_direction, "minimize")) {
    exposure_preference <- -exposure_preference
  }
  ranks <- .pareto_rank_2d(
    criterion_preference,
    exposure_preference
  )
  depth_score <- .pareto_depth_score(ranks)

  pareto_rank <- rep(NA_integer_, n)
  pareto_depth_score <- rep(NA_real_, n)
  relative_priority <- rep(NA_real_, n)
  pareto_rank[included] <- ranks
  pareto_depth_score[included] <- depth_score
  relative_priority[included] <- depth_score
  front_sizes <- tabulate(ranks, nbins = max(ranks))

  table <- data.frame(
    cell = seq_len(n),
    reference_weight = reference_weight,
    pareto_rank = pareto_rank,
    pareto_depth_score = pareto_depth_score,
    relative_priority = relative_priority,
    included = included,
    stringsAsFactors = FALSE
  )
  if (!identical(criterion_data$field, "reference_weight")) {
    table[[criterion_data$field]] <- criterion_data$values
  }
  table[[exposure_data$field]] <- exposure_data$values
  table <- table[, c(
    "cell",
    "reference_weight",
    if (!identical(criterion_data$field, "reference_weight")) {
      criterion_data$field
    },
    exposure_data$field,
    "pareto_rank",
    "pareto_depth_score",
    "relative_priority",
    "included"
  ), drop = FALSE]

  out <- list(
    call = match.call(),
    exposure = exposure_data$field,
    exposure_label = exposure_data$label,
    criterion = criterion_data$field,
    criterion_label = criterion_data$label,
    criterion_direction = criterion_direction,
    exposure_direction = exposure_direction,
    scope = scope,
    positive_only = positive_only,
    table = table,
    front_sizes = data.frame(
      pareto_rank = seq_along(front_sizes),
      n = front_sizes
    ),
    rasters = .priority_rasters(
      x,
      pareto_rank = pareto_rank,
      pareto_depth_score = pareto_depth_score
    )
  )
  class(out) <- "climniche_priority"
  out
}

#' @rdname climniche_priority
#' @param object A `climniche_priority` object.
#' @param ... Unused.
#' @return `summary()` returns a `summary.climniche_priority` object with the
#'   fitted settings and Pareto diagnostics.
#' @export
summary.climniche_priority <- function(object, ...) {
  if (!inherits(object, "climniche_priority")) {
    stop("object must be a climniche_priority object.", call. = FALSE)
  }
  included <- object$table$included
  exposure <- object$table[[object$exposure]][included]
  exposure_direction <- object$exposure_direction %||% "maximize"
  if (identical(exposure_direction, "minimize")) {
    exposure <- -exposure
  }
  criterion <- object$table[[object$criterion]][included]
  if (identical(object$criterion_direction, "minimize")) {
    criterion <- -criterion
  }
  correlation <- if (length(unique(exposure)) > 1L &&
                     length(unique(criterion)) > 1L) {
    stats::cor(exposure, criterion, method = "spearman")
  } else {
    NA_real_
  }
  first_front_n <- object$front_sizes$n[1L]
  diagnostics <- data.frame(
    ranked_cells = sum(included),
    pareto_fronts = nrow(object$front_sizes),
    first_front_cells = first_front_n,
    first_front_fraction = first_front_n / sum(included),
    objective_rank_correlation = correlation,
    unique_exposure_values = length(unique(exposure)),
    unique_criterion_values = length(unique(criterion)),
    stringsAsFactors = FALSE
  )
  out <- list(
    settings = data.frame(
      exposure = object$exposure_label,
      exposure_direction = exposure_direction,
      criterion = object$criterion_label,
      criterion_direction = object$criterion_direction,
      scope = object$scope,
      positive_only = object$positive_only,
      stringsAsFactors = FALSE
    ),
    diagnostics = diagnostics,
    front_sizes = object$front_sizes
  )
  class(out) <- "summary.climniche_priority"
  out
}

#' @export
print.summary.climniche_priority <- function(x, ...) {
  cat("Ecological screening of climate exposure\n\n")
  cat("Exposure:", x$settings$exposure, "\n")
  cat("Exposure direction:", x$settings$exposure_direction, "\n")
  cat("Second criterion:", x$settings$criterion, "\n")
  cat("Criterion direction:", x$settings$criterion_direction, "\n")
  cat("Scope:", x$settings$scope, "\n")
  cat("Positive exposure only:", x$settings$positive_only, "\n\n")
  cat("Ranked cells:", x$diagnostics$ranked_cells, "\n")
  cat("Pareto fronts:", x$diagnostics$pareto_fronts, "\n")
  cat(
    "First-front cells:", x$diagnostics$first_front_cells,
    sprintf("(%.1f%%)\n", 100 * x$diagnostics$first_front_fraction)
  )
  cat(
    "Spearman correlation between objectives:",
    format(x$diagnostics$objective_rank_correlation, digits = 3),
    "\n"
  )
  invisible(x)
}

#' @export
print.climniche_priority <- function(x, ...) {
  included <- x$table$included
  first_front_n <- x$front_sizes$n[1L]
  cat("climniche ecological screening of climate exposure\n\n")
  cat("Exposure:", x$exposure_label, "\n")
  cat("Exposure direction:", x$exposure_direction %||% "maximize", "\n")
  cat("Second criterion:", x$criterion_label,
      paste0(" (", x$criterion_direction, ")\n"))
  cat("Scope:", x$scope, "\n")
  cat("Positive exposure only:", x$positive_only, "\n")
  cat("Ranked cells:", sum(included), "\n")
  cat("Pareto fronts:", nrow(x$front_sizes), "\n")
  cat(
    "First-front cells:", first_front_n,
    sprintf("(%.1f%%)\n", 100 * first_front_n / sum(included))
  )
  invisible(x)
}

.priority_exposure_values <- function(x, exposure, positive_only) {
  key <- .metric_key(exposure)
  choices <- c(
    "niche_distance_change",
    "niche_boundary_exceedance",
    "climate_reconfiguration",
    "climate_change_amount"
  )
  if (!key %in% choices) {
    stop("Unsupported exposure quantity.", call. = FALSE)
  }

  values <- .fit_metric(x, key)
  if (!is.numeric(values) || length(values) != nrow(x$current)) {
    stop("x does not contain the requested exposure quantity.", call. = FALSE)
  }
  values <- as.numeric(values)

  labels <- c(
    climate_change_amount = "Climatic Displacement",
    niche_distance_change = "Niche Distance Shift",
    climate_reconfiguration = "Climatic Reconfiguration",
    niche_boundary_exceedance = "Niche Boundary Exceedance"
  )
  list(values = values, field = key, label = unname(labels[[key]]))
}

.priority_criterion_values <- function(x, criterion, criterion_name = NULL) {
  if (is.null(criterion)) {
    values <- .fit_reference_weights(x)
    field <- "reference_weight"
    label <- criterion_name %||% "Reference weight"
    return(list(values = values, field = field, label = label))
  }

  values <- .priority_values_from_input(x, criterion)
  field <- "decision_criterion"
  label <- criterion_name %||% "Decision criterion"
  list(values = values, field = field, label = label)
}

.priority_values_from_input <- function(x, criterion) {
  n <- nrow(x$current)
  if (is.numeric(criterion)) {
    values <- criterion
    if (length(values) == n) {
      return(as.numeric(.align_reference_values(
        values,
        rownames(x$current),
        arg = "criterion"
      )))
    }
    if (!is.null(x$raster_complete) &&
        length(values) == length(x$raster_complete)) {
      return(as.numeric(values[x$raster_complete]))
    }
    stop("criterion must match the evaluated rows or raster cells.",
         call. = FALSE)
  }

  template <- .fit_spatial_template(x)
  if (is.null(template) || is.null(x$raster_complete)) {
    stop("A raster criterion requires a spatial climniche fit.", call. = FALSE)
  }
  if (.is_raster(criterion)) {
    .need_raster()
    if (!.is_raster(template) || raster::nlayers(criterion) != 1L ||
        !raster::compareRaster(template, criterion, stopiffalse = FALSE)) {
      stop("criterion raster must match the fitted raster geometry.",
           call. = FALSE)
    }
    values <- raster::getValues(criterion)
  } else if (.is_spatraster(criterion)) {
    .need_terra()
    if (!.is_spatraster(template) || terra::nlyr(criterion) != 1L ||
        !terra::compareGeom(template, criterion, stopOnError = FALSE)) {
      stop("criterion raster must match the fitted raster geometry.",
           call. = FALSE)
    }
    values <- terra::values(criterion)[, 1L]
  } else {
    stop("criterion must be numeric or a matching one-layer raster.",
         call. = FALSE)
  }
  as.numeric(values[x$raster_complete])
}

.pareto_rank_2d <- function(first, second) {
  if (!is.numeric(first) || !is.numeric(second) ||
      length(first) != length(second) || !length(first) ||
      any(!is.finite(first)) || any(!is.finite(second))) {
    stop("Pareto criteria must be finite numeric vectors of equal length.",
         call. = FALSE)
  }

  second_levels <- sort(unique(second), decreasing = TRUE)
  second_index <- match(second, second_levels)
  order_index <- order(first, second, decreasing = TRUE)
  tree <- integer(length(second_levels))

  query_tree <- function(index) {
    best <- 0L
    while (index > 0L) {
      if (tree[index] > best) best <- tree[index]
      index <- index - bitwAnd(index, -index)
    }
    best
  }
  update_tree <- function(index, value) {
    n_tree <- length(tree)
    while (index <= n_tree) {
      if (value > tree[index]) tree[index] <<- value
      index <- index + bitwAnd(index, -index)
    }
  }

  out <- integer(length(first))
  position <- 1L
  n <- length(order_index)
  while (position <= n) {
    end <- position
    while (end < n &&
           first[order_index[end + 1L]] == first[order_index[position]] &&
           second[order_index[end + 1L]] == second[order_index[position]]) {
      end <- end + 1L
    }
    rows <- order_index[position:end]
    rank <- query_tree(second_index[rows[1L]]) + 1L
    out[rows] <- rank
    update_tree(second_index[rows[1L]], rank)
    position <- end + 1L
  }
  out
}

.pareto_depth_score <- function(ranks) {
  maximum <- max(ranks)
  if (maximum == 1L) {
    return(rep(1, length(ranks)))
  }
  (maximum - ranks) / (maximum - 1)
}

.priority_rasters <- function(x, pareto_rank, pareto_depth_score) {
  rasters <- list(
    pareto_rank = .priority_value_raster(x, pareto_rank),
    pareto_depth_score = .priority_value_raster(x, pareto_depth_score),
    relative_priority = .priority_value_raster(x, pareto_depth_score)
  )
  if (all(vapply(rasters, is.null, logical(1)))) return(NULL)
  for (name in names(rasters)) names(rasters[[name]]) <- name
  rasters
}

.priority_value_raster <- function(x, values) {
  template <- .fit_spatial_template(x)
  complete <- x$raster_complete
  if (is.null(template) || is.null(complete)) return(NULL)
  if (length(values) != sum(complete)) {
    stop("Screening values do not match the fitted spatial cells.",
         call. = FALSE)
  }
  if (.is_raster(template)) {
    return(.values_to_raster(template, values, complete))
  }
  if (.is_spatraster(template)) {
    return(.values_to_spatraster(template, values, complete))
  }
  NULL
}

#' Plot ecological screening of climate exposure
#'
#' @param x A `climniche_priority` object.
#' @param type `"plane"`, `"map"`, or `"both"`.
#' @param map_value Map `"pareto_depth_score"`, the original `"pareto_rank"`,
#'   or the legacy alias `"relative_priority"`.
#' @param max_points Target maximum number of cells drawn in the decision plane.
#'   Pareto rank 1 is always retained.
#' @param seed Random seed used when the decision plane is subsampled.
#' @param title Optional overall title for a combined figure.
#' @param extent Optional `c(xmin, xmax, ymin, ymax)` map extent.
#' @param degree_labels Longitude-latitude label style.
#' @param study_region Optional study-region boundary accepted by
#'   [plot_climniche_map()].
#' @param legend_position Legend position for the map.
#' @param ... Additional arguments passed to `plot_climniche_priority()` by the
#'   S3 `plot()` method.
#'
#' @return A ggplot object, a patchwork object, or a named list of plots when
#'   `patchwork` is unavailable.
#'
#' @examples
#' sim <- simulate_climniche(n = 400, p = 6, seed = 21)
#' fit <- fit_climniche(
#'   sim[["current"]],
#'   sim[["future_away"]],
#'   occupied = sim[["occupied"]],
#'   sensitivity = sim[["sensitivity"]]
#' )
#' priority <- climniche_priority(fit)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_climniche_priority(priority, type = "plane")
#' }
#' @export
plot_climniche_priority <- function(
    x,
    type = c("plane", "map", "both"),
    map_value = c("pareto_depth_score", "pareto_rank", "relative_priority"),
    max_points = 6000L,
    seed = 1L,
    title = NULL,
    extent = NULL,
    degree_labels = c("auto", "none", "hemisphere"),
    study_region = NULL,
    legend_position = "bottom") {
  .need_ggplot2()
  if (!inherits(x, "climniche_priority")) {
    stop("x must be a climniche_priority object.", call. = FALSE)
  }
  type <- match.arg(type)
  map_value <- match.arg(map_value)
  degree_labels <- match.arg(degree_labels)
  max_points <- .check_positive_integer(max_points, "max_points")
  seed <- .check_finite_scalar(seed, "seed")
  if (seed != floor(seed)) {
    stop("seed must be an integer.", call. = FALSE)
  }

  plane <- .plot_priority_plane(x, max_points = max_points, seed = seed)
  if (identical(type, "plane")) return(plane)

  map <- .plot_priority_map(
    x,
    map_value = map_value,
    extent = extent,
    degree_labels = degree_labels,
    study_region = study_region,
    legend_position = legend_position
  )
  if (identical(type, "map")) return(map)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(list(plane = plane, map = map))
  }
  plane <- plane +
    ggplot2::labs(title = "(a) Exposure and ecological criterion") +
    ggplot2::theme(legend.position = "none")
  map <- map + ggplot2::labs(
    title = paste0("(b) ", .priority_map_title(x, map_value))
  )
  combined <- plane + map + patchwork::plot_layout(widths = c(0.85, 1.35))
  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(
      title = title
    )
  }
  combined
}

#' @rdname plot_climniche_priority
#' @export
plot.climniche_priority <- function(x, ...) {
  plot_climniche_priority(x, ...)
}

.priority_colours <- function(exposure_direction = "maximize") {
  if (identical(exposure_direction, "minimize")) {
    return(c("#f2f5f4", "#c9dfda", "#82b8ae", "#3a8078", "#104a46"))
  }
  c("#f6f2e8", "#f2ce85", "#e49a48", "#c75b4d", "#7d2431")
}

.plot_priority_plane <- function(x, max_points, seed) {
  criterion_field <- x$criterion
  exposure_field <- x$exposure
  exposure_direction <- x$exposure_direction %||% "maximize"
  depth_field <- if ("pareto_depth_score" %in% names(x$table)) {
    "pareto_depth_score"
  } else {
    "relative_priority"
  }
  data <- x$table[x$table$included, c(
    criterion_field, exposure_field, "pareto_rank", depth_field
  ), drop = FALSE]
  names(data) <- c(
    "criterion_value", "exposure_value", "pareto_rank", "pareto_depth_score"
  )
  front <- data[data$pareto_rank == 1L, , drop = FALSE]

  if (nrow(data) > max_points) {
    keep_front <- which(data$pareto_rank == 1L)
    available <- setdiff(seq_len(nrow(data)), keep_front)
    n_sample <- max(0L, max_points - length(keep_front))
    set.seed(as.integer(seed))
    sampled <- if (n_sample < length(available)) {
      sample(available, n_sample)
    } else {
      available
    }
    data <- data[sort(c(keep_front, sampled)), , drop = FALSE]
  }

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = criterion_value,
      y = exposure_value,
      colour = pareto_depth_score
    )
  ) +
    ggplot2::geom_point(size = 0.65, alpha = 0.48) +
    ggplot2::geom_point(
      data = front,
      mapping = ggplot2::aes(x = criterion_value, y = exposure_value),
      shape = 21,
      fill = NA,
      colour = "black",
      size = 1.15,
      stroke = 0.25,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_gradientn(
      colours = .priority_colours(exposure_direction),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = "Exposure and ecological criterion",
      subtitle = .priority_plane_subtitle(x),
      x = x$criterion_label,
      y = x$exposure_label,
      colour = "Pareto depth score"
    ) +
    .climniche_theme(base_size = 8.5) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = grid::unit(12, "mm")
    )
}

.plot_priority_map <- function(x, map_value, extent, degree_labels,
                               study_region, legend_position) {
  raster_name <- map_value
  if (identical(map_value, "pareto_depth_score") &&
      is.null(x$rasters[[map_value]]) &&
      !is.null(x$rasters[["relative_priority"]])) {
    raster_name <- "relative_priority"
  }
  if (is.null(x$rasters) || is.null(x$rasters[[raster_name]])) {
    stop("x does not contain spatial Pareto layers.", call. = FALSE)
  }
  raster <- x$rasters[[raster_name]]
  data <- .spatial_df(raster, value = "value")
  limits <- .map_extent(data, extent = extent)
  cell_size <- .spatial_res(raster)
  axis_spec <- .map_axis_spec(raster, degree_labels = degree_labels)
  region <- .study_region_df(study_region, raster)

  if (map_value %in% c("pareto_depth_score", "relative_priority")) {
    colours <- .priority_colours(x$exposure_direction %||% "maximize")
    scale_limits <- c(0, 1)
    legend_title <- "Pareto depth score"
  } else {
    colours <- rev(.priority_colours(
      x$exposure_direction %||% "maximize"
    ))
    scale_limits <- range(data$value, na.rm = TRUE)
    if (diff(scale_limits) == 0) {
      scale_limits <- scale_limits + c(-0.5, 0.5)
    }
    legend_title <- "Pareto rank"
  }

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = x, y = y, fill = value)
  ) +
    ggplot2::geom_tile(width = cell_size[1L], height = cell_size[2L]) +
    ggplot2::coord_equal(
      xlim = limits$xlim,
      ylim = limits$ylim,
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(labels = axis_spec$x) +
    ggplot2::scale_y_continuous(labels = axis_spec$y) +
    ggplot2::scale_fill_gradientn(
      colours = colours,
      limits = scale_limits,
      na.value = NA,
      oob = .squish_values
    ) +
    ggplot2::labs(
      title = .priority_map_title(x, map_value),
      x = axis_spec$xlab,
      y = axis_spec$ylab,
      fill = legend_title
    ) +
    .map_theme() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.direction = "horizontal",
      legend.key.width = grid::unit(12, "mm")
    )
  if (!is.null(region) && nrow(region) > 0L) {
    plot <- plot + ggplot2::geom_path(
      data = region,
      ggplot2::aes(x = x, y = y, group = group),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.3,
      lineend = "round"
    )
  }
  plot
}

.priority_plane_subtitle <- function(x) {
  exposure_direction <- x$exposure_direction %||% "maximize"
  criterion_direction <- x$criterion_direction %||% "maximize"
  exposure_text <- if (identical(exposure_direction, "maximize")) {
    "higher exposure"
  } else {
    "lower exposure"
  }
  criterion_text <- if (identical(criterion_direction, "maximize")) {
    "higher criterion values"
  } else {
    "lower criterion values"
  }
  paste0(
    "Pareto rank 1 is outlined; ", exposure_text, " and ", criterion_text,
    " are preferred"
  )
}

.priority_map_title <- function(x, map_value) {
  switch(
    map_value,
    pareto_depth_score = "Spatial Pareto depth",
    relative_priority = "Spatial Pareto depth",
    pareto_rank = "Spatial Pareto rank"
  )
}
