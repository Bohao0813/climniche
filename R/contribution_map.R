#' Summarise dominant climatic contributions by cell
#'
#' Identifies the climate variable with the largest absolute contribution to
#' future minus current squared niche distance at each evaluated cell.
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` retains cells with positive current reference
#'   weight; `"all"` retains every evaluated cell.
#'
#' @return A `climniche_contribution` object. `table` contains the dominant
#'   variable, signed contribution and dominance share for each cell. `summary`
#'   contains mean contribution shares and signed contributions by variable,
#'   together with the fraction of non-zero analysis weight for which each
#'   variable is uniquely dominant. `squared_niche_distance_change` is the
#'   row sum of the signed variable contributions. The legacy field
#'   `niche_potential_change` is retained as an exact alias. Spatial fits also
#'   return raster layers.
#'
#' @details
#' Let \eqn{V_{ij}} be the contribution returned by
#' [variable_contribution()] for cell \eqn{i} and variable \eqn{j}. The dominant
#' variable has the largest \eqn{|V_{ij}|}. Its dominance share is
#' \deqn{H_i = \frac{\max_j |V_{ij}|}{\sum_j |V_{ij}|}.}
#' The share is undefined when every contribution is zero. Equal largest
#' absolute contributions are reported as `"Tied"` rather than being assigned
#' according to column order.
#'
#' Contributions sum exactly to future minus current squared niche distance.
#' This total is also
#' \eqn{(r_{1i} - r_{0i})(r_{1i} + r_{0i})}. The terms therefore attribute the
#' squared-distance change underlying Niche Distance Shift. For a non-diagonal
#' metric matrix, each variable term includes its part of the cross-variable
#' terms in the chosen climatic basis. They are not SDM variable importance or
#' causal effects.
#' Current-scope summaries use the fitted reference weights; all-scope
#' summaries give each evaluated cell equal weight.
#'
#' @examples
#' sim <- simulate_climniche(n = 300, p = 6, seed = 31)
#' fit <- fit_climniche(
#'   sim[["current"]],
#'   sim[["future_away"]],
#'   occupied = sim[["occupied"]],
#'   sensitivity = sim[["sensitivity"]]
#' )
#' contribution <- climniche_dominant_contribution(fit)
#' summary(contribution)
#' head(contribution[["table"]])
#' @export
climniche_dominant_contribution <- function(
    x,
    scope = c("current", "all")) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  values <- .as_numeric_matrix(
    x$variable_contribution,
    "x$variable_contribution"
  )
  if (nrow(values) != nrow(x$current)) {
    stop("x contains incompatible variable contributions.", call. = FALSE)
  }
  variables <- colnames(values)
  if (is.null(variables) || any(!nzchar(variables)) || anyDuplicated(variables)) {
    variables <- paste0("climate_", seq_len(ncol(values)))
  }

  reference_weight <- .fit_reference_weights(x)
  included <- rep(TRUE, nrow(values))
  if (identical(scope, "current")) included <- reference_weight > 0

  absolute <- abs(values)
  total_absolute <- rowSums(absolute)
  maximum_absolute <- apply(absolute, 1L, max)
  tolerance <- sqrt(.Machine$double.eps) * pmax(1, maximum_absolute)
  ties <- rowSums(abs(absolute - maximum_absolute) <= tolerance)
  dominant_index <- max.col(absolute, ties.method = "first")
  defined <- included & is.finite(total_absolute) & total_absolute > 0
  unique_dominant <- defined & ties == 1L
  tied <- defined & ties > 1L

  dominant_variable <- rep(NA_character_, nrow(values))
  dominant_variable[unique_dominant] <- variables[
    dominant_index[unique_dominant]
  ]
  dominant_variable[tied] <- "Tied"
  dominant_id <- rep(NA_integer_, nrow(values))
  dominant_id[unique_dominant] <- dominant_index[unique_dominant]
  dominant_id[tied] <- length(variables) + 1L

  dominant_contribution <- rep(NA_real_, nrow(values))
  dominant_contribution[unique_dominant] <- values[cbind(
    which(unique_dominant),
    dominant_index[unique_dominant]
  )]
  dominant_share <- rep(NA_real_, nrow(values))
  dominant_share[defined] <- maximum_absolute[defined] /
    total_absolute[defined]
  squared_niche_distance_change <- rowSums(values)

  cell <- if (!is.null(x$cell_index) &&
              length(x$cell_index) == nrow(values)) {
    x$cell_index
  } else {
    seq_len(nrow(values))
  }
  table <- data.frame(
    cell = cell,
    reference_weight = reference_weight,
    included = included,
    dominant_variable = dominant_variable,
    dominant_variable_id = dominant_id,
    dominant_contribution = dominant_contribution,
    dominant_share = dominant_share,
    tied = tied,
    total_absolute_contribution = total_absolute,
    squared_niche_distance_change = squared_niche_distance_change,
    niche_potential_change = squared_niche_distance_change,
    stringsAsFactors = FALSE
  )

  analysis_weight <- if (identical(scope, "current")) {
    reference_weight
  } else {
    as.numeric(included)
  }
  share <- matrix(NA_real_, nrow = nrow(values), ncol = ncol(values))
  valid_share <- total_absolute > 0
  share[valid_share, ] <- absolute[valid_share, , drop = FALSE] /
    total_absolute[valid_share]
  analysis_rows <- included & is.finite(analysis_weight) & analysis_weight > 0
  if (any(defined)) {
    mean_share <- .weighted_col_means(
      share[defined, , drop = FALSE],
      analysis_weight[defined]
    )
    dominant_fraction <- vapply(seq_len(ncol(values)), function(j) {
      .weighted_mean_vector(
        as.numeric(dominant_id[defined] == j),
        analysis_weight[defined]
      )
    }, numeric(1))
    mean_dominant_share <- vapply(seq_len(ncol(values)), function(j) {
      selected <- defined & dominant_id == j
      if (!any(selected)) return(NA_real_)
      .weighted_mean_vector(
        dominant_share[selected],
        analysis_weight[selected]
      )
    }, numeric(1))
  } else {
    mean_share <- dominant_fraction <- mean_dominant_share <-
      rep(NA_real_, ncol(values))
  }
  if (any(analysis_rows)) {
    mean_signed <- .weighted_col_means(
      values[analysis_rows, , drop = FALSE],
      analysis_weight[analysis_rows]
    )
    positive_fraction <- vapply(seq_len(ncol(values)), function(j) {
      .weighted_mean_vector(
        as.numeric(values[analysis_rows, j] > 0),
        analysis_weight[analysis_rows]
      )
    }, numeric(1))
  } else {
    mean_signed <- positive_fraction <- rep(NA_real_, ncol(values))
  }
  summary <- data.frame(
    variable = variables,
    mean_absolute_share = as.numeric(mean_share),
    mean_signed_contribution = as.numeric(mean_signed),
    positive_contribution_fraction = positive_fraction,
    dominant_weight_fraction = dominant_fraction,
    mean_dominant_share = mean_dominant_share,
    stringsAsFactors = FALSE
  )
  summary <- summary[order(
    summary$mean_absolute_share,
    decreasing = TRUE
  ), , drop = FALSE]
  rownames(summary) <- NULL

  raster_values <- list(
    dominant_variable = dominant_id,
    dominant_share = dominant_share,
    dominant_contribution = dominant_contribution,
    squared_niche_distance_change = squared_niche_distance_change,
    niche_potential_change = squared_niche_distance_change
  )
  raster_values <- lapply(raster_values, function(value) {
    value[!included] <- NA_real_
    value
  })
  rasters <- lapply(raster_values, function(value) {
    .contribution_value_raster(x, value)
  })
  rasters <- rasters[!vapply(rasters, is.null, logical(1))]
  for (name in names(rasters)) names(rasters[[name]]) <- name

  out <- list(
    call = match.call(),
    scope = scope,
    variables = variables,
    table = table,
    summary = summary,
    lookup = data.frame(
      id = seq_len(length(variables) + 1L),
      variable = c(variables, "Tied"),
      stringsAsFactors = FALSE
    ),
    rasters = if (length(rasters)) rasters else NULL
  )
  class(out) <- "climniche_contribution"
  out
}

#' @rdname climniche_dominant_contribution
#' @param object A `climniche_contribution` object.
#' @param ... Additional arguments passed to methods.
#' @export
summary.climniche_contribution <- function(object, ...) {
  if (!inherits(object, "climniche_contribution")) {
    stop("object must be a climniche_contribution object.", call. = FALSE)
  }
  object$summary
}

#' @export
print.climniche_contribution <- function(x, ...) {
  included <- x$table$included
  defined <- included & !is.na(x$table$dominant_variable)
  cat("Dominant climatic contributions\n\n")
  cat("Scope:", x$scope, "\n")
  cat("Analysed cells:", sum(included), "\n")
  cat("Cells with non-zero contributions:", sum(defined), "\n")
  cat("Variables:", length(x$variables), "\n")
  invisible(x)
}

.contribution_value_raster <- function(x, values) {
  template <- .fit_spatial_template(x)
  complete <- x$raster_complete
  if (is.null(template) || is.null(complete)) return(NULL)
  if (length(values) != sum(complete)) {
    stop("Contribution values do not match the fitted spatial cells.",
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

#' Map dominant climatic contributions
#'
#' @param x A `climniche_contribution` object or fitted `climniche_fit` object.
#' @param type Plot the dominant `"variable"`, its `"share"`, or `"both"`.
#' @param scope Scope used when `x` is a fitted object.
#' @param variable_labels Optional named vector replacing climate variable
#'   names.
#' @param legend_variables Optional character vector of fitted variable names
#'   to retain in the legend, including variables that are not dominant in any
#'   mapped cell. The default shows only observed categories.
#' @param colours Optional colours for the climate variables. A final colour
#'   may be supplied for tied contributions; otherwise ties are grey.
#' @param title Optional overall title for a combined figure.
#' @param extent Optional `c(xmin, xmax, ymin, ymax)` map extent.
#' @param degree_labels Longitude-latitude label style.
#' @param study_region Optional study-region boundary accepted by
#'   [plot_climniche_map()].
#' @param legend_position Legend position.
#' @param ... Additional arguments passed to
#'   `plot_climniche_dominant_contribution()` by the S3 `plot()` method.
#'
#' @return A ggplot object, a patchwork object, or a named list when patchwork
#'   is unavailable.
#'
#' @examples
#' sim <- simulate_climniche(n = 300, p = 6, seed = 32)
#' fit <- fit_climniche(
#'   sim[["current"]],
#'   sim[["future_away"]],
#'   occupied = sim[["occupied"]],
#'   sensitivity = sim[["sensitivity"]]
#' )
#' contribution <- climniche_dominant_contribution(fit)
#' @export
plot_climniche_dominant_contribution <- function(
    x,
    type = c("variable", "share", "both"),
    scope = c("current", "all"),
    variable_labels = NULL,
    colours = NULL,
    title = NULL,
    extent = NULL,
    degree_labels = c("auto", "none", "hemisphere"),
    study_region = NULL,
    legend_position = "bottom",
    legend_variables = NULL) {
  .need_ggplot2()
  type <- match.arg(type)
  scope <- match.arg(scope)
  degree_labels <- match.arg(degree_labels)
  if (inherits(x, "climniche_fit")) {
    x <- climniche_dominant_contribution(x, scope = scope)
  }
  if (!inherits(x, "climniche_contribution")) {
    stop("x must be a climniche_contribution or climniche_fit object.",
         call. = FALSE)
  }
  if (is.null(x$rasters)) {
    stop("x does not contain spatial contribution layers.", call. = FALSE)
  }

  labels <- x$lookup$variable
  if (!is.null(variable_labels)) {
    if (is.null(names(variable_labels)) || any(!nzchar(names(variable_labels)))) {
      stop("variable_labels must be a named vector.", call. = FALSE)
    }
    matched <- match(labels, names(variable_labels))
    replace <- !is.na(matched)
    labels[replace] <- unname(variable_labels[matched[replace]])
  }
  legend_breaks <- ggplot2::waiver()
  drop_legend_levels <- TRUE
  if (!is.null(legend_variables)) {
    legend_variables <- unique(as.character(legend_variables))
    if (!length(legend_variables) || anyNA(legend_variables) ||
        any(!nzchar(legend_variables))) {
      stop("legend_variables must contain fitted variable names.",
           call. = FALSE)
    }
    matched <- match(legend_variables, x$lookup$variable)
    if (anyNA(matched)) {
      stop("legend_variables contains names not found in x.", call. = FALSE)
    }
    legend_breaks <- labels[matched]
    drop_legend_levels <- FALSE
  }
  if (is.null(colours)) {
    colours <- .contribution_colours(nrow(x$lookup))
  }
  if (length(colours) == length(x$variables)) {
    colours <- c(colours, "#777777")
  }
  if (length(colours) < nrow(x$lookup)) {
    stop("colours must contain one value per climate variable.",
         call. = FALSE)
  }
  colours <- colours[seq_len(nrow(x$lookup))]
  names(colours) <- labels

  variable_plot <- .plot_dominant_variable_map(
    x,
    labels = labels,
    colours = colours,
    legend_breaks = legend_breaks,
    drop_legend_levels = drop_legend_levels,
    extent = extent,
    degree_labels = degree_labels,
    study_region = study_region,
    legend_position = legend_position
  )
  if (identical(type, "variable")) return(variable_plot)
  share_plot <- .plot_dominant_share_map(
    x,
    extent = extent,
    degree_labels = degree_labels,
    study_region = study_region,
    legend_position = legend_position
  )
  if (identical(type, "share")) return(share_plot)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(list(variable = variable_plot, share = share_plot))
  }
  variable_plot <- variable_plot + ggplot2::labs(
    title = "(a) Dominant climatic contribution"
  )
  share_plot <- share_plot + ggplot2::labs(title = "(b) Dominance share")
  combined <- variable_plot + share_plot +
    patchwork::plot_layout(widths = c(1.15, 1))
  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(title = title)
  }
  combined
}

#' @rdname plot_climniche_dominant_contribution
#' @export
plot.climniche_contribution <- function(x, ...) {
  plot_climniche_dominant_contribution(x, ...)
}

.contribution_colours <- function(n) {
  base <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#6A3D9A", "#8C6D31", "#000000"
  )
  if (n <= 1L) return("#777777")
  variable_n <- n - 1L
  variable_colours <- if (variable_n <= length(base)) {
    base[seq_len(variable_n)]
  } else {
    grDevices::hcl.colors(variable_n, palette = "Dark 3")
  }
  c(variable_colours, "#777777")
}

.contribution_map_data <- function(raster, value, extent, degree_labels,
                                   study_region) {
  data <- .spatial_df(raster, value = value)
  list(
    data = data,
    limits = .map_extent(data, extent = extent),
    cell_size = .spatial_res(raster),
    axis = .map_axis_spec(raster, degree_labels = degree_labels),
    region = .study_region_df(study_region, raster)
  )
}

.add_contribution_region <- function(plot, region) {
  if (is.null(region) || !nrow(region)) return(plot)
  plot + ggplot2::geom_path(
    data = region,
    ggplot2::aes(x = x, y = y, group = group),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.3,
    lineend = "round"
  )
}

.plot_dominant_variable_map <- function(x, labels, colours, legend_breaks,
                                        drop_legend_levels, extent,
                                        degree_labels, study_region,
                                        legend_position) {
  map <- .contribution_map_data(
    x$rasters$dominant_variable,
    value = "variable_id",
    extent = extent,
    degree_labels = degree_labels,
    study_region = study_region
  )
  map$data$variable <- factor(
    map$data$variable_id,
    levels = x$lookup$id,
    labels = labels
  )
  plot <- ggplot2::ggplot(
    map$data,
    ggplot2::aes(x = x, y = y, fill = variable)
  ) +
    ggplot2::geom_tile(
      width = map$cell_size[1L],
      height = map$cell_size[2L],
      show.legend = TRUE
    ) +
    ggplot2::coord_equal(
      xlim = map$limits$xlim,
      ylim = map$limits$ylim,
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(labels = map$axis$x) +
    ggplot2::scale_y_continuous(labels = map$axis$y) +
    ggplot2::scale_fill_manual(
      values = colours,
      breaks = legend_breaks,
      drop = drop_legend_levels,
      na.value = NA
    ) +
    ggplot2::labs(
      title = "Dominant climatic contribution",
      x = map$axis$xlab,
      y = map$axis$ylab,
      fill = "Climate variable"
    ) +
    .map_theme() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_text(colour = "black")
    )
  .add_contribution_region(plot, map$region)
}

.plot_dominant_share_map <- function(x, extent, degree_labels, study_region,
                                     legend_position) {
  map <- .contribution_map_data(
    x$rasters$dominant_share,
    value = "dominant_share",
    extent = extent,
    degree_labels = degree_labels,
    study_region = study_region
  )
  plot <- ggplot2::ggplot(
    map$data,
    ggplot2::aes(x = x, y = y, fill = dominant_share)
  ) +
    ggplot2::geom_tile(
      width = map$cell_size[1L],
      height = map$cell_size[2L]
    ) +
    ggplot2::coord_equal(
      xlim = map$limits$xlim,
      ylim = map$limits$ylim,
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(labels = map$axis$x) +
    ggplot2::scale_y_continuous(labels = map$axis$y) +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7f7f4", "#bfd5cf", "#70a89f", "#2d7069", "#123f3c"),
      limits = c(0, 1),
      na.value = NA,
      oob = .squish_values
    ) +
    ggplot2::labs(
      title = "Dominance share",
      x = map$axis$xlab,
      y = map$axis$ylab,
      fill = "Dominance share"
    ) +
    .map_theme() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.direction = "horizontal",
      legend.key.width = grid::unit(12, "mm"),
      legend.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_text(colour = "black")
    )
  .add_contribution_region(plot, map$region)
}
