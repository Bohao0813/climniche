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

.dominant_exposure_bins <- function(tab, bins = 45L) {
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
    stringsAsFactors = FALSE
  )
  if (!nrow(dat)) {
    return(data.frame(
      x_mid = numeric(), y_mid = numeric(), class = factor(),
      class_count = integer(), total_count = integer(),
      proportion_in_bin = numeric(), cell_weight = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  counts <- stats::aggregate(
    list(class_count = rep(1L, nrow(dat))),
    by = list(x_id = dat$x_id, y_id = dat$y_id, class = dat$class),
    FUN = sum
  )
  counts <- counts[order(counts$x_id, counts$y_id, -counts$class_count,
                         counts$class), , drop = FALSE]
  bin_key <- paste(counts$x_id, counts$y_id, sep = ":")
  total <- stats::aggregate(
    list(total_count = counts$class_count),
    by = list(x_id = counts$x_id, y_id = counts$y_id),
    FUN = sum
  )
  total_key <- paste(total$x_id, total$y_id, sep = ":")
  out <- counts[!duplicated(bin_key), , drop = FALSE]
  out_key <- paste(out$x_id, out$y_id, sep = ":")
  out$total_count <- total$total_count[match(out_key, total_key)]
  out$proportion_in_bin <- out$class_count / out$total_count
  out$x_mid <- (x_breaks[out$x_id] + x_breaks[out$x_id + 1L]) / 2
  out$y_mid <- (y_breaks[out$y_id] + y_breaks[out$y_id + 1L]) / 2
  out$x_width <- diff(x_breaks)[out$x_id]
  out$y_height <- diff(y_breaks)[out$y_id]
  out$cell_weight <- sqrt(out$total_count / max(out$total_count, na.rm = TRUE))
  out$class <- factor(out$class, levels = names(.class_colours()))
  out[, c("x_mid", "y_mid", "class", "class_count", "total_count",
          "proportion_in_bin", "cell_weight", "x_width", "y_height")]
}

.showcase_metric_distribution <- function(tab) {
  metrics <- c(
    climate_change_amount = "Climate change\namount",
    niche_distance_change = "Niche distance\nchange",
    composition_change = "Composition\nchange",
    outside_niche_exceedance = "Niche boundary\nexceedance"
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

#' Build data for the climniche showcase figure
#'
#' @param x A fitted climniche object.
#' @param scope `"current"` for current occurrence, range or thresholded SDM
#'   cells; `"all"` for all evaluated cells.
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
  tab$class <- factor(as.character(tab$class), levels = names(.class_colours()))

  plane <- tab
  if (nrow(plane) > max_points) {
    set.seed(seed)
    plane <- plane[sample(seq_len(nrow(plane)), max_points), , drop = FALSE]
  }

  class_counts <- as.data.frame(table(tab$class), stringsAsFactors = FALSE)
  names(class_counts) <- c("class", "count")
  class_counts$proportion <- class_counts$count / nrow(tab)
  class_counts <- class_counts[class_counts$count > 0, , drop = FALSE]
  class_counts$class <- factor(class_counts$class,
                               levels = rev(as.character(class_counts$class)))
  class_counts$label <- ifelse(
    class_counts$proportion > 0 & class_counts$proportion < 0.01,
    "<1%",
    paste0(round(100 * class_counts$proportion), "%")
  )

  idx <- if (scope == "current") x$occupied else seq_len(nrow(x$current))
  vals <- colMeans(x$variable_contribution[idx, , drop = FALSE],
                   na.rm = TRUE)
  variables <- data.frame(
    variable = names(vals),
    mean_contribution = as.numeric(vals),
    abs_mean_contribution = abs(as.numeric(vals)),
    direction = ifelse(vals >= 0,
                       "less similar to current niche",
                       "more similar to current niche"),
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
    b_potential <- as.numeric(stats::quantile(
      x$psi_current[x$occupied],
      probs = q,
      names = FALSE,
      na.rm = TRUE,
      type = 8
    ))
    b_radius <- sqrt(pmax(0, b_potential))
    exceed <- pmax(0, x$niche_radius_future[idx] - b_radius)
    data.frame(
      boundary_quantile = q,
      boundary_distance = b_radius,
      prop_exceeded = mean(exceed > 0, na.rm = TRUE),
      mean_exceedance = mean(exceed, na.rm = TRUE)
    )
  }))

  out <- list(
    plane = plane,
    plane_bins = .dominant_exposure_bins(tab, bins = plane_bins),
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

#' Plot the climniche showcase figure
#'
#' @param x A fitted climniche object or data returned by
#'   `climniche_showcase_data()`.
#' @param scope `"current"` for current occurrence, range or thresholded SDM
#'   cells; `"all"` for all evaluated cells.
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
    "little climate niche change" = "little climate\nniche change",
    "closer to current niche" = "closer to\ncurrent niche",
    "farther from current niche" = "farther from\ncurrent niche",
    "outside current niche boundary" = "outside current\nniche boundary",
    "changed composition, similar distance" = "changed composition,\nsimilar distance"
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
      x = "Climate change amount",
      y = "Niche distance change",
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
    ggplot2::labs(
      title = "(b) Exposure classes",
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
      labels = c("TRUE" = "less similar", "FALSE" = "more similar")
    ) +
    ggplot2::labs(
      title = "(c) Variable contributions",
      x = "Mean contribution to niche distance change",
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
