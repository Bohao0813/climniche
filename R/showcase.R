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

.weighted_quantity_bins <- function(tab, x, y, bins = 45L, weights = NULL) {
  bins <- max(8L, as.integer(bins)[1])
  x_values <- tab[[x]]
  y_values <- tab[[y]]
  x_breaks <- .showcase_breaks(x_values, bins)
  y_breaks <- .showcase_breaks(y_values, bins)
  if (is.null(weights)) {
    weights <- if ("occupied_weight" %in% names(tab)) {
      tab$occupied_weight
    } else {
      rep(1, nrow(tab))
    }
  }
  ok <- is.finite(x_values) & is.finite(y_values) &
    is.finite(weights) & weights > 0
  if (!any(ok)) {
    return(data.frame(
      x_mid = numeric(), y_mid = numeric(), total_count = integer(),
      total_weight = numeric(), x_width = numeric(), y_height = numeric()
    ))
  }
  dat <- data.frame(
    x_id = findInterval(x_values[ok], x_breaks, all.inside = TRUE),
    y_id = findInterval(y_values[ok], y_breaks, all.inside = TRUE),
    weight = weights[ok]
  )
  out <- stats::aggregate(
    list(total_count = rep(1L, nrow(dat)), total_weight = dat$weight),
    by = list(x_id = dat$x_id, y_id = dat$y_id),
    FUN = sum
  )
  out$x_mid <- (x_breaks[out$x_id] + x_breaks[out$x_id + 1L]) / 2
  out$y_mid <- (y_breaks[out$y_id] + y_breaks[out$y_id + 1L]) / 2
  out$x_width <- diff(x_breaks)[out$x_id]
  out$y_height <- diff(y_breaks)[out$y_id]
  out[, c("x_mid", "y_mid", "total_count", "total_weight",
          "x_width", "y_height")]
}

.showcase_metric_distribution <- function(tab, weights) {
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
      weight = weights,
      stringsAsFactors = FALSE
    )
  }))
  out <- out[is.finite(out$value) & is.finite(out$weight) & out$weight > 0,
             , drop = FALSE]
  out$metric <- factor(out$metric, levels = unname(metrics))
  out
}

.showcase_metric_histograms <- function(metric_data, bins = 30L) {
  bins <- max(10L, as.integer(bins)[1])
  metric_levels <- levels(metric_data$metric)
  histograms <- list()
  zero_mass <- data.frame(
    metric = factor(character(), levels = metric_levels),
    proportion = numeric(),
    label = character()
  )

  for (metric_name in metric_levels) {
    dat <- metric_data[metric_data$metric == metric_name, , drop = FALSE]
    total_weight <- sum(dat$weight)
    if (!nrow(dat) || total_weight <= 0) {
      next
    }
    is_boundary <- identical(metric_name, "Niche Boundary\nExceedance")
    zero <- if (is_boundary) dat$value <= sqrt(.Machine$double.eps) else FALSE
    if (is_boundary) {
      proportion <- sum(dat$weight[zero]) / total_weight
      zero_mass <- rbind(
        zero_mass,
        data.frame(
          metric = factor(metric_name, levels = metric_levels),
          proportion = proportion,
          label = paste0("No exceedance: ", round(100 * proportion), "%")
        )
      )
      dat <- dat[!zero, , drop = FALSE]
    }
    if (!nrow(dat)) {
      next
    }
    breaks <- .showcase_breaks(dat$value, bins)
    bin <- findInterval(dat$value, breaks, all.inside = TRUE)
    weighted_count <- tapply(dat$weight, bin, sum)
    ids <- as.integer(names(weighted_count))
    histograms[[length(histograms) + 1L]] <- data.frame(
      metric = factor(metric_name, levels = metric_levels),
      x = (breaks[ids] + breaks[ids + 1L]) / 2,
      width = diff(breaks)[ids],
      proportion = as.numeric(weighted_count) / total_weight,
      stringsAsFactors = FALSE
    )
  }

  histogram_data <- if (length(histograms)) {
    do.call(rbind, histograms)
  } else {
    data.frame(
      metric = factor(character(), levels = metric_levels),
      x = numeric(), width = numeric(), proportion = numeric()
    )
  }

  list(
    histograms = histogram_data,
    zero_mass = zero_mass
  )
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
                                    plane_bins = 35L,
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
  class_counts$class <- factor(class_counts$class,
                               levels = rev(as.character(class_counts$class)))
  class_counts$label <- ifelse(
    class_counts$proportion == 0,
    "0%",
    ifelse(class_counts$proportion < 0.01,
    "<1%",
    paste0(round(100 * class_counts$proportion), "%"))
  )
  descriptors <- .descriptor_summary(tab, tab_weight)
  descriptor_levels <- c(
    "Away from realised niche centre",
    "Limited Niche Distance Shift",
    "Toward realised niche centre",
    "Beyond empirical niche boundary",
    "Within empirical niche boundary"
  )
  descriptors$level <- factor(descriptors$level,
                              levels = descriptor_levels)
  descriptors$label <- ifelse(
    descriptors$proportion == 0,
    "0%",
    ifelse(descriptors$proportion < 0.01, "<1%",
           paste0(round(100 * descriptors$proportion), "%"))
  )

  idx <- if (scope == "current") x$occupied else seq_len(nrow(x$current))
  reference_weights <- .fit_reference_weights(x)
  weights <- if (scope == "current") {
    reference_weights[idx]
  } else {
    rep(1, length(idx))
  }
  contribution <- x$variable_contribution[idx, , drop = FALSE]
  vals <- .weighted_col_means(contribution, weights)
  absolute_contribution <- abs(contribution)
  contribution_total <- rowSums(absolute_contribution)
  contribution_share <- absolute_contribution
  positive_total <- contribution_total > 0
  contribution_share[positive_total, ] <-
    contribution_share[positive_total, , drop = FALSE] /
    contribution_total[positive_total]
  contribution_share[!positive_total, ] <- 0
  mean_absolute_share <- .weighted_col_means(contribution_share, weights)
  variables <- data.frame(
    variable = names(vals),
    mean_contribution = as.numeric(vals),
    abs_mean_contribution = abs(as.numeric(vals)),
    mean_absolute_share = as.numeric(mean_absolute_share),
    direction = ifelse(vals >= 0,
                       "positive niche potential contribution",
                       "negative niche potential contribution"),
    stringsAsFactors = FALSE
  )
  variables <- variables[order(variables$mean_absolute_share,
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

  metric_data <- .showcase_metric_distribution(tab, tab_weight)
  metric_histograms <- .showcase_metric_histograms(metric_data)

  out <- list(
    plane = plane,
    plane_bins = .dominant_exposure_bins(tab, bins = plane_bins,
                                         weights = tab_weight),
    reconfiguration_bins = .weighted_quantity_bins(
      tab,
      x = "climate_reconfiguration",
      y = "niche_boundary_exceedance",
      bins = plane_bins,
      weights = tab_weight
    ),
    classes = class_counts,
    descriptors = descriptors,
    variables = variables,
    boundary = boundary,
    metrics = metric_data,
    metric_histograms = metric_histograms$histograms,
    metric_zero_mass = metric_histograms$zero_mass,
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

#' Build data for the climniche summary figure
#'
#' `climniche_summary_figure_data()` is a neutral alias for
#' [climniche_showcase_data()]. The older function name is retained for
#' compatibility.
#'
#' @param x A fitted climniche object.
#' @param ... Arguments passed to [climniche_showcase_data()].
#'
#' @return A list of data frames used by [plot_climniche_summary_figure()].
#' @export
climniche_summary_figure_data <- function(x, ...) {
  climniche_showcase_data(x, ...)
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
                                    plane_bins = 35L,
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

  base_size <- 7.0
  tile_width <- if (nrow(x$plane_bins)) x$plane_bins$x_width[1] else 1
  tile_height <- if (nrow(x$plane_bins)) x$plane_bins$y_height[1] else 1
  p_plane <- ggplot2::ggplot(
    x$plane_bins,
    ggplot2::aes(x = x_mid, y = y_mid, fill = total_weight)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25,
                        linetype = 2, colour = "grey45") +
    ggplot2::geom_tile(width = tile_width, height = tile_height,
                       colour = NA) +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7fbfa", "#b7d7d1", "#4f8f86", "#1f5f58"),
      trans = "sqrt",
      name = "Summed\nreference weight"
    ) +
    ggplot2::labs(
      title = "(a) Climatic Displacement and Niche Distance Shift",
      x = "Climatic Displacement",
      y = "Niche Distance Shift"
    ) +
    .climniche_theme(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        title.position = "top",
        barwidth = grid::unit(2.8, "mm"),
        barheight = grid::unit(23, "mm")
      )
    )

  secondary_width <- if (nrow(x$reconfiguration_bins)) {
    x$reconfiguration_bins$x_width[1]
  } else {
    1
  }
  secondary_height <- if (nrow(x$reconfiguration_bins)) {
    x$reconfiguration_bins$y_height[1]
  } else {
    1
  }
  p_reconfiguration <- ggplot2::ggplot(
    x$reconfiguration_bins,
    ggplot2::aes(x = x_mid, y = y_mid, fill = total_weight)
  ) +
    ggplot2::geom_tile(
      width = secondary_width,
      height = secondary_height,
      colour = NA
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#fffaf0", "#f2ce91", "#d99238", "#8a3f20"),
      trans = "sqrt",
      name = "Summed\nreference weight"
    ) +
    ggplot2::labs(
      title = "(b) Climatic Reconfiguration and Niche Boundary Exceedance",
      x = "Climatic Reconfiguration",
      y = "Niche Boundary Exceedance"
    ) +
    .climniche_theme(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        title.position = "top",
        barwidth = grid::unit(2.8, "mm"),
        barheight = grid::unit(23, "mm")
      )
    )

  contribution_direction <- x$variables$mean_contribution >= 0
  mixed_direction <- length(unique(contribution_direction)) > 1L
  p_vars <- ggplot2::ggplot(
    x$variables,
    ggplot2::aes(x = mean_absolute_share, y = variable)
  )
  if (mixed_direction) {
    p_vars <- p_vars +
      ggplot2::geom_col(
        ggplot2::aes(fill = mean_contribution > 0),
        width = 0.66, colour = "grey25", linewidth = 0.15
      ) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#c65d57", "FALSE" = "#4c78a8"),
        labels = c("TRUE" = "Positive", "FALSE" = "Negative")
      )
  } else {
    p_vars <- p_vars + ggplot2::geom_col(
      width = 0.66, fill = "#6f9fbc",
      colour = "grey25", linewidth = 0.15
    )
  }
  direction_subtitle <- if (mixed_direction) {
    "Colour indicates the direction of the fitted variable contribution"
  } else {
    NULL
  }
  p_vars <- p_vars +
    ggplot2::scale_x_continuous(
      labels = function(z) paste0(round(100 * z), "%")
    ) +
    ggplot2::labs(
      title = "(c) Variable Contributions",
      subtitle = direction_subtitle,
      x = "Mean absolute contribution",
      y = NULL
    ) +
    .climniche_theme(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (mixed_direction) "bottom" else "none",
      legend.direction = "horizontal",
      legend.key.width = grid::unit(3.6, "mm"),
      legend.key.height = grid::unit(3.0, "mm")
    )
  if (mixed_direction) {
    p_vars <- p_vars + ggplot2::labs(fill = "Fitted contribution")
  }

  p_metrics <- ggplot2::ggplot(
    x$metric_histograms,
    ggplot2::aes(x = x, y = proportion)
  ) +
    ggplot2::geom_col(ggplot2::aes(width = width),
                      fill = "#6f9fbc", colour = "white",
                      linewidth = 0.10) +
    ggplot2::geom_vline(
      data = data.frame(
        metric = factor("Niche Distance\nShift",
                        levels = levels(x$metric_histograms$metric)),
        xintercept = 0
      ),
      ggplot2::aes(xintercept = xintercept),
      linewidth = 0.22, linetype = 2, colour = "grey45"
    ) +
    ggplot2::geom_text(
      data = x$metric_zero_mass,
      ggplot2::aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = 1.02, vjust = 1.15,
      size = 1.9, colour = "black"
    ) +
    ggplot2::facet_wrap(~metric, scales = "free_x", nrow = 2) +
    ggplot2::scale_y_continuous(
      labels = function(z) paste0(round(100 * z), "%")
    ) +
    ggplot2::labs(
      title = "(d) Metric Distributions",
      x = NULL,
      y = "Weighted percentage"
    ) +
    .climniche_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "none")

  plots <- list(
    exposure = p_plane,
    reconfiguration = p_reconfiguration,
    variables = p_vars,
    metrics = p_metrics
  )
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_vars_layout <- p_vars
    if ("free" %in% getNamespaceExports("patchwork")) {
      free_args <- names(formals(patchwork::free))
      p_vars_layout <- if (all(c("type", "side") %in% free_args)) {
        patchwork::free(p_vars, type = "space", side = "l")
      } else {
        patchwork::free(p_vars)
      }
    }
    out <- (p_plane | p_reconfiguration) / (p_vars_layout | p_metrics) +
      patchwork::plot_layout(widths = c(1, 1.18))
    if (!is.null(title)) {
      out <- out + patchwork::plot_annotation(title = title)
    }
    return(out)
  }
  plots
}

#' Plot the climniche summary figure
#'
#' `plot_climniche_summary_figure()` is a neutral alias for
#' [plot_climniche_showcase()]. The older function name is retained for
#' compatibility.
#'
#' @param x A fitted climniche object or data returned by
#'   [climniche_showcase_data()].
#' @param ... Arguments passed to [plot_climniche_showcase()].
#'
#' @return A patchwork object when `patchwork` is installed, otherwise a named
#'   list of ggplot objects.
#' @export
plot_climniche_summary_figure <- function(x, ...) {
  plot_climniche_showcase(x, ...)
}
