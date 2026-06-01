.showcase_breaks <- function(x, bins) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(seq(0, 1, length.out = bins + 1L))
  }
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) <= 0) {
    pad <- if (rng[1] == 0) 0.5 else abs(rng[1]) * 0.05
    rng <- rng + c(-pad, pad)
  } else {
    pad <- diff(rng) * 0.001
    rng <- rng + c(-pad, pad)
  }
  seq(rng[1], rng[2], length.out = bins + 1L)
}

.dominant_exposure_bins <- function(tab, bins = 45L, weights = NULL) {
  bins <- max(8L, as.integer(bins)[1])
  x_breaks <- .showcase_breaks(tab$climate_change_amount, bins)
  y_breaks <- .showcase_breaks(tab$niche_distance_change, bins)
  x_id <- findInterval(tab$climate_change_amount, x_breaks,
                       all.inside = TRUE)
  y_id <- findInterval(tab$niche_distance_change, y_breaks,
                       all.inside = TRUE)
  ok <- is.finite(tab$climate_change_amount) &
    is.finite(tab$niche_distance_change) &
    !is.na(tab$class)
  dat <- data.frame(
    x_id = x_id[ok],
    y_id = y_id[ok],
    class = as.character(tab$class[ok]),
    occupied_weight = if (!is.null(weights)) {
      weights[ok]
    } else if ("occupied_weight" %in% names(tab)) {
      tab$occupied_weight[ok]
    } else {
      rep(1, sum(ok))
    },
    stringsAsFactors = FALSE
  )
  if (!nrow(dat)) {
    return(data.frame(
      x_mid = numeric(), y_mid = numeric(), class = factor(),
      class_count = integer(), total_count = integer(),
      class_weight = numeric(), total_weight = numeric(),
      proportion_in_bin = numeric(), cell_weight = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  counts <- stats::aggregate(
    list(class_count = rep(1L, nrow(dat))),
    by = list(x_id = dat$x_id, y_id = dat$y_id, class = dat$class),
    FUN = sum
  )
  class_weight <- stats::aggregate(
    list(class_weight = dat$occupied_weight),
    by = list(x_id = dat$x_id, y_id = dat$y_id, class = dat$class),
    FUN = sum
  )
  counts <- merge(counts, class_weight,
                  by = c("x_id", "y_id", "class"), sort = FALSE)
  counts <- counts[order(counts$x_id, counts$y_id, -counts$class_weight,
                         counts$class), , drop = FALSE]
  bin_key <- paste(counts$x_id, counts$y_id, sep = ":")
  total <- stats::aggregate(
    list(total_count = counts$class_count),
    by = list(x_id = counts$x_id, y_id = counts$y_id),
    FUN = sum
  )
  total_weight <- stats::aggregate(
    list(total_weight = counts$class_weight),
    by = list(x_id = counts$x_id, y_id = counts$y_id),
    FUN = sum
  )
  total_key <- paste(total$x_id, total$y_id, sep = ":")
  out <- counts[!duplicated(bin_key), , drop = FALSE]
  out_key <- paste(out$x_id, out$y_id, sep = ":")
  out$total_count <- total$total_count[match(out_key, total_key)]
  out$total_weight <- total_weight$total_weight[match(out_key, total_key)]
  out$proportion_in_bin <- out$class_weight / out$total_weight
  out$x_mid <- (x_breaks[out$x_id] + x_breaks[out$x_id + 1L]) / 2
  out$y_mid <- (y_breaks[out$y_id] + y_breaks[out$y_id + 1L]) / 2
  out$x_width <- diff(x_breaks)[out$x_id]
  out$y_height <- diff(y_breaks)[out$y_id]
  out$cell_weight <- sqrt(out$total_weight / max(out$total_weight, na.rm = TRUE))
  out$class <- factor(out$class, levels = names(.class_colours()))
  out[, c("x_mid", "y_mid", "class", "class_count", "total_count",
          "class_weight", "total_weight", "proportion_in_bin",
          "cell_weight", "x_width", "y_height")]
}

.showcase_metric_distribution <- function(tab) {
  metrics <- c(
    climate_change_amount = "Climatic\nDisplacement",
    niche_distance_change = "Niche Distance\nShift",
    climate_reconfiguration = "Climatic\nReconfiguration",
    niche_boundary_exceedance = "Niche Boundary\nExceedance"
  )
  out <- do.call(rbind, lapply(names(metrics), function(nm) {
    data.frame(
      metric = unname(metrics[nm]),
      value = tab[[nm]],
      stringsAsFactors = FALSE
    )
  }))
  out <- out[is.finite(out$value), , drop = FALSE]
  out$metric <- factor(out$metric, levels = unname(metrics))
  out
}

#' Build data for the climniche summary figure
#'
#' @param x A fitted climniche object.
#' @param scope `"current"` for current reference cells; `"all"` for all
#'   evaluated cells.
#' @param max_points Maximum number of cells to keep for the exposure plane.
#' @param seed Random seed used when subsampling cells.
#' @param plane_bins Number of fixed bins used to summarize the exposure plane.
#' @param boundary_probs Boundary quantiles used for the sensitivity curve.
#' @param top_variables Number of variables to show.
#'
#' @return A list of data frames used by `plot_climniche_showcase()`.
#' @export
climniche_showcase_data <- function(x, scope = c("current", "all"),
                                    max_points = 6000L, seed = 1L,
                                    plane_bins = 45L,
                                    boundary_probs = seq(0.50, 0.99, 0.01),
                                    top_variables = 6L) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  tab <- climniche_table(x, scope = scope)
  tab$class <- .normalise_class(tab$class)
  tab_weight <- if (scope == "current") tab$occupied_weight else rep(1, nrow(tab))

  plane <- tab
  if (nrow(plane) > max_points) {
    set.seed(seed)
    plane <- plane[sample(seq_len(nrow(plane)), max_points), , drop = FALSE]
  }

  class_counts <- as.data.frame(table(tab$class), stringsAsFactors = FALSE)
  names(class_counts) <- c("class", "count")
  class_weight <- tapply(tab_weight, tab$class, sum)
  class_counts$weight <- as.numeric(class_weight[as.character(class_counts$class)])
  class_counts$weight[is.na(class_counts$weight)] <- 0
  class_counts$proportion <- class_counts$weight / sum(tab_weight)
  class_counts <- class_counts[class_counts$count > 0, , drop = FALSE]
  class_counts$class <- factor(class_counts$class,
                               levels = rev(as.character(class_counts$class)))
  class_counts$label <- ifelse(
    class_counts$proportion > 0 & class_counts$proportion < 0.01,
    "<1%",
    paste0(round(100 * class_counts$proportion), "%")
  )

  idx <- if (scope == "current") x$occupied else seq_len(nrow(x$current))
  reference_weights <- .fit_reference_weights(x)
  weights <- if (scope == "current") {
    reference_weights[idx]
  } else {
    rep(1, length(idx))
  }
  vals <- .weighted_col_means(x$variable_contribution[idx, , drop = FALSE],
                              weights)
  variables <- data.frame(
    variable = names(vals),
    mean_contribution = as.numeric(vals),
    abs_mean_contribution = abs(as.numeric(vals)),
    direction = ifelse(vals >= 0,
                       "positive niche potential contribution",
                       "negative niche potential contribution"),
    stringsAsFactors = FALSE
  )
  variables <- variables[order(variables$abs_mean_contribution,
                               decreasing = TRUE), , drop = FALSE]
  variables <- utils::head(variables, top_variables)
  variables$variable <- factor(variables$variable,
                               levels = rev(variables$variable))

  boundary_probs <- sort(unique(boundary_probs))
  boundary_probs <- boundary_probs[is.finite(boundary_probs) &
                                     boundary_probs > 0 & boundary_probs < 1]
  boundary <- do.call(rbind, lapply(boundary_probs, function(q) {
    b_potential <- as.numeric(.weighted_quantile(
      x$psi_current,
      weights = reference_weights,
      probs = q,
      names = FALSE
    ))
    b_radius <- sqrt(pmax(0, b_potential))
    exceed <- pmax(0, x$niche_radius_future[idx] - b_radius)
    data.frame(
      boundary_quantile = q,
      boundary_distance = b_radius,
      prop_exceeded = .weighted_prop(
        exceed > x$classification_settings$boundary_exceedance_tolerance,
        weights
      ),
      mean_exceedance = .weighted_mean_vector(exceed, weights)
    )
  }))

  out <- list(
    plane = plane,
    plane_bins = .dominant_exposure_bins(tab, bins = plane_bins,
                                         weights = tab_weight),
    classes = class_counts,
    variables = variables,
    boundary = boundary,
    metrics = .showcase_metric_distribution(tab),
    settings = data.frame(
      scope = scope,
      n_cells = nrow(tab),
      sampled_cells = nrow(plane),
      plane_bins = plane_bins,
      fitted_boundary_quantile = x$boundary_quantile,
      fitted_boundary_distance = x$boundary_radius
    )
  )
  class(out) <- "climniche_showcase_data"
  out
}

#' Plot the climniche summary figure
#'
#' @param x A fitted climniche object or data returned by
#'   `climniche_showcase_data()`.
#' @param scope `"current"` for current reference cells; `"all"` for all
#'   evaluated cells.
#' @param max_points Maximum number of cells to draw in the exposure plane.
#' @param seed Random seed used when subsampling cells.
#' @param plane_bins Number of fixed bins used to summarize the exposure plane.
#' @param boundary_probs Boundary quantiles used for the sensitivity curve.
#' @param top_variables Number of variables to show.
#' @param variable_labels Optional named vector replacing variable labels.
#' @param title Optional overall title when `patchwork` is installed.
#'
#' @return A patchwork object when `patchwork` is installed, otherwise a named
#'   list of ggplot objects.
#' @export
plot_climniche_showcase <- function(x, scope = c("current", "all"),
                                    max_points = 6000L, seed = 1L,
                                    plane_bins = 45L,
                                    boundary_probs = seq(0.50, 0.99, 0.01),
                                    top_variables = 6L,
                                    variable_labels = NULL,
                                    title = NULL) {
  .need_ggplot2()
  if (!inherits(x, "climniche_showcase_data")) {
    x <- climniche_showcase_data(
      x,
      scope = scope,
      max_points = max_points,
      seed = seed,
      plane_bins = plane_bins,
      boundary_probs = boundary_probs,
      top_variables = top_variables
    )
  }
  if (!is.null(variable_labels)) {
    idx <- match(as.character(x$variables$variable), names(variable_labels))
    replace <- !is.na(idx)
    levels(x$variables$variable)[match(
      as.character(x$variables$variable)[replace],
      levels(x$variables$variable)
    )] <- unname(variable_labels[idx[replace]])
  }

  class_cols <- .class_colours()
  class_labs <- c(
    "Limited climate niche change" = "Limited climate\nniche change",
    "Closer to current niche" = "Closer to\ncurrent niche",
    "Farther from current niche" = "Farther from\ncurrent niche",
    "Outside current niche boundary" = "Outside current\nniche boundary",
    "Climatic Reconfiguration with limited Niche Distance Shift" =
      "Climatic Reconfiguration\nwith limited Niche Distance Shift"
  )
  tile_width <- if (nrow(x$plane_bins)) x$plane_bins$x_width[1] else 1
  tile_height <- if (nrow(x$plane_bins)) x$plane_bins$y_height[1] else 1
  p_plane <- ggplot2::ggplot(
    x$plane_bins,
    ggplot2::aes(x = x_mid, y = y_mid, fill = class,
                 alpha = cell_weight)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25,
                        linetype = 2, colour = "grey45") +
    ggplot2::geom_tile(width = tile_width, height = tile_height,
                       colour = NA) +
    ggplot2::scale_fill_manual(
      values = class_cols,
      labels = class_labs,
      drop = TRUE,
      guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::scale_alpha(range = c(0.25, 0.95), guide = "none") +
    ggplot2::labs(
      title = "(a) Binned exposure plane",
      subtitle = "Squares are fixed bins; colours follow panel (b)",
      x = "Climatic Displacement",
      y = "Niche Distance Shift",
      fill = "Dominant exposure class"
    ) +
    .climniche_theme(base_size = 8.2) +
    ggplot2::theme(
      legend.position = "none"
    )

  p_classes <- ggplot2::ggplot(
    x$classes,
    ggplot2::aes(x = proportion, y = class, fill = class)
  ) +
    ggplot2::geom_col(width = 0.66, colour = "grey25", linewidth = 0.15) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.12,
                       size = 2.3) +
    ggplot2::scale_x_continuous(
      labels = function(z) paste0(round(100 * z), "%"),
      limits = c(0, max(x$classes$proportion) * 1.18),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::scale_fill_manual(values = class_cols, drop = TRUE) +
    ggplot2::scale_y_discrete(labels = class_labs) +
    ggplot2::labs(
      title = "(b) Derived exposure classes",
      x = "Proportion of analysed cells",
      y = NULL
    ) +
    .climniche_theme(base_size = 8.2) +
    ggplot2::theme(legend.position = "none")

  p_vars <- ggplot2::ggplot(
    x$variables,
    ggplot2::aes(x = mean_contribution, y = variable,
                 fill = mean_contribution > 0)
  ) +
    ggplot2::geom_col(width = 0.66, colour = "grey25", linewidth = 0.15) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25,
                        colour = "grey45") +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#c65d57", "FALSE" = "#4c78a8"),
      labels = c("TRUE" = "positive", "FALSE" = "negative")
    ) +
    ggplot2::labs(
      title = "(c) Variable contributions",
      x = "Mean contribution to niche potential change",
      y = NULL,
      fill = NULL
    ) +
    .climniche_theme(base_size = 8.2) +
    ggplot2::theme(legend.position = "none")

  p_metrics <- ggplot2::ggplot(
    x$metrics,
    ggplot2::aes(x = value)
  ) +
    ggplot2::geom_histogram(bins = 32, fill = "#6f9fbc",
                            colour = "white", linewidth = 0.10) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.22,
                        linetype = 2, colour = "grey45") +
    ggplot2::facet_wrap(~metric, scales = "free", nrow = 2) +
    ggplot2::labs(
      title = "(d) Metric distributions",
      x = NULL,
      y = "Number of cells"
    ) +
    .climniche_theme(base_size = 8.2) +
    ggplot2::theme(legend.position = "none")

  plots <- list(
    exposure = p_plane,
    classes = p_classes,
    variables = p_vars,
    metrics = p_metrics
  )
  if (requireNamespace("patchwork", quietly = TRUE)) {
    out <- (p_plane | p_classes) / (p_vars | p_metrics)
    if (!is.null(title)) {
      out <- out + patchwork::plot_annotation(title = title)
    }
    return(out)
  }
  plots
}
