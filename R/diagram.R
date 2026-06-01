.weighted_niche_coordinates <- function(x) {
  ev <- eigen(x$A, symmetric = TRUE)
  keep <- ev$values > sqrt(.Machine$double.eps)
  if (!any(keep)) {
    stop("The niche metric has no positive dimensions.", call. = FALSE)
  }
  transform <- ev$vectors[, keep, drop = FALSE] %*%
    diag(sqrt(ev$values[keep]), nrow = sum(keep))
  current <- sweep(x$current, 2L, x$center, "-") %*% transform
  future <- sweep(x$future, 2L, x$center, "-") %*% transform
  list(current = current, future = future, transform = transform)
}

.niche_diagram_rotation <- function(coords, weights) {
  idx <- .positive_reference_indices(weights)
  ref <- coords[idx, , drop = FALSE]
  weights <- weights[idx]
  if (ncol(ref) == 1L) {
    return(matrix(c(1, 0), nrow = 1))
  }
  S <- crossprod(ref * sqrt(weights / sum(weights)))
  ev <- eigen(S, symmetric = TRUE)
  rot <- ev$vectors[, seq_len(min(2L, ncol(ev$vectors))), drop = FALSE]
  if (ncol(rot) == 1L) {
    rot <- cbind(rot, 0)
  }
  rot
}

.hull_df <- function(dat, x = "axis1", y = "axis2") {
  ok <- is.finite(dat[[x]]) & is.finite(dat[[y]])
  dat <- dat[ok, , drop = FALSE]
  if (nrow(dat) < 3L) {
    return(dat[0, , drop = FALSE])
  }
  dat[grDevices::chull(dat[[x]], dat[[y]]), , drop = FALSE]
}

.circle_df <- function(radius, n = 240L) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    axis1 = radius * cos(theta),
    axis2 = radius * sin(theta)
  )
}

#' Build data for a niche climate exposure diagram
#'
#' @param x A fitted climniche object.
#' @param scope `"current"` for current reference cells; `"all"` for all
#'   evaluated cells.
#' @param max_arrows Maximum number of current to future arrows to keep.
#' @param seed Random seed used when subsampling arrows.
#'
#' @return A list of data frames used by `plot_climniche_diagram()`.
#' @export
climniche_diagram_data <- function(x, scope = c("current", "all"),
                                   max_arrows = 350L, seed = 1L) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  coords <- .weighted_niche_coordinates(x)
  weights_all <- if (!is.null(x$occupied_weight)) {
    x$occupied_weight
  } else {
    replace(rep(0, nrow(x$current)), x$occupied, 1)
  }
  rot <- .niche_diagram_rotation(coords$current, weights_all)
  current2 <- coords$current %*% rot
  future2 <- coords$future %*% rot
  colnames(current2) <- colnames(future2) <- c("axis1", "axis2")

  tab <- climniche_table(x, scope = scope)
  idx <- tab$cell
  cells <- data.frame(
    cell = idx,
    current_axis1 = current2[idx, 1],
    current_axis2 = current2[idx, 2],
    future_axis1 = future2[idx, 1],
    future_axis2 = future2[idx, 2],
    climate_change_amount = tab$climate_change_amount,
    niche_distance_change = tab$niche_distance_change,
    climate_reconfiguration = tab$climate_reconfiguration,
    composition_change = tab$climate_reconfiguration,
    niche_boundary_exceedance = tab$niche_boundary_exceedance,
    outside_niche_exceedance = tab$niche_boundary_exceedance,
    occupied_weight = tab$occupied_weight,
    class = tab$class,
    mixed_variable_response = tab$mixed_variable_response,
    stringsAsFactors = FALSE
  )
  cells$class <- .normalise_class(cells$class)
  class_summary <- do.call(rbind, lapply(levels(cells$class), function(cls) {
    rows <- cells$class == cls
    rows[is.na(rows)] <- FALSE
    if (!any(rows)) {
      return(NULL)
    }
    w <- if (scope == "current") cells$occupied_weight[rows] else rep(1, sum(rows))
    data.frame(
      class = cls,
      current_axis1 = .weighted_mean_vector(cells$current_axis1[rows], w),
      current_axis2 = .weighted_mean_vector(cells$current_axis2[rows], w),
      future_axis1 = .weighted_mean_vector(cells$future_axis1[rows], w),
      future_axis2 = .weighted_mean_vector(cells$future_axis2[rows], w),
      climate_change_amount =
        .weighted_mean_vector(cells$climate_change_amount[rows], w),
      niche_distance_change =
        .weighted_mean_vector(cells$niche_distance_change[rows], w),
      climate_reconfiguration =
        .weighted_mean_vector(cells$climate_reconfiguration[rows], w),
      composition_change =
        .weighted_mean_vector(cells$climate_reconfiguration[rows], w),
      niche_boundary_exceedance =
        .weighted_mean_vector(cells$niche_boundary_exceedance[rows], w),
      outside_niche_exceedance =
        .weighted_mean_vector(cells$niche_boundary_exceedance[rows], w),
      stringsAsFactors = FALSE
    )
  }))
  class_counts <- as.data.frame(table(cells$class), stringsAsFactors = FALSE)
  names(class_counts) <- c("class", "count")
  class_summary <- merge(class_summary, class_counts, by = "class",
                         all.x = TRUE, sort = FALSE)
  class_weights <- tapply(
    if (scope == "current") cells$occupied_weight else rep(1, nrow(cells)),
    cells$class,
    sum
  )
  class_summary$proportion <- as.numeric(
    class_weights[as.character(class_summary$class)] /
      sum(class_weights, na.rm = TRUE)
  )
  class_summary$class <- .normalise_class(class_summary$class)

  arrows <- cells
  if (nrow(arrows) > max_arrows) {
    set.seed(seed)
    arrows <- arrows[sample(seq_len(nrow(arrows)), max_arrows), , drop = FALSE]
  }

  reference <- data.frame(axis1 = current2[, 1], axis2 = current2[, 2])
  current_ref <- data.frame(axis1 = current2[x$occupied, 1],
                            axis2 = current2[x$occupied, 2])
  future_ref <- data.frame(axis1 = future2[idx, 1], axis2 = future2[idx, 2])

  basis <- diag(ncol(x$current)) %*% coords$transform %*% rot
  if (is.null(colnames(x$current))) {
    var_names <- paste0("var", seq_len(nrow(basis)))
  } else {
    var_names <- colnames(x$current)
  }
  arrow_radius <- stats::quantile(
    sqrt(rowSums(rbind(current2, future2)^2)),
    probs = 0.90,
    na.rm = TRUE,
    names = FALSE
  )
  basis_len <- sqrt(rowSums(basis^2))
  scale <- if (all(basis_len == 0)) 1 else 0.85 * arrow_radius / max(basis_len)
  variables <- data.frame(
    variable = var_names,
    axis1 = basis[, 1] * scale,
    axis2 = basis[, 2] * scale,
    stringsAsFactors = FALSE
  )

  settings <- x$classification_settings
  if (is.null(settings)) {
    settings <- list()
  }

  out <- list(
    cells = cells,
    arrows = arrows,
    class_summary = class_summary,
    reference_hull = .hull_df(reference),
    current_hull = .hull_df(current_ref),
    future_hull = .hull_df(future_ref),
    niche_boundary = .circle_df(x$boundary_radius),
    variables = variables,
    center = data.frame(axis1 = 0, axis2 = 0),
    settings = data.frame(
      scope = scope,
      boundary_quantile = x$boundary_quantile,
      boundary_radius = x$boundary_radius,
      tolerance = settings$tolerance %||% NA_real_,
      stable_climate_change =
        settings$stable_climate_change %||% NA_real_,
      stable_reconfiguration =
        settings$stable_reconfiguration %||% NA_real_,
      boundary_exceedance_tolerance =
        settings$boundary_exceedance_tolerance %||% NA_real_,
      n_cells = nrow(cells),
      n_arrows = nrow(arrows)
    )
  )
  class(out) <- "climniche_diagram_data"
  out
}

#' Plot a niche climate exposure diagram
#'
#' @noRd
.plot_climniche_metric_diagram <- function(x, show_reference = FALSE,
                                           show_hulls = TRUE,
                                           boundary_shape = c("hull", "circle", "none"),
                                           show_boundary_label = TRUE,
                                           show_center = TRUE,
                                           title = NULL,
                                           max_points = 6000L,
                                           seed = 1L) {
  boundary_shape <- match.arg(boundary_shape)
  cell_cols <- .class_colours()

  cells <- x$cells
  if (nrow(cells) > max_points) {
    set.seed(seed)
    cells <- cells[sample(seq_len(nrow(cells)), max_points), , drop = FALSE]
  }

  p_space <- ggplot2::ggplot()
  if (show_reference && nrow(x$reference_hull) > 2L) {
    p_space <- p_space + ggplot2::geom_polygon(
      data = x$reference_hull,
      ggplot2::aes(x = axis1, y = axis2),
      fill = "#f5f5f5", colour = "#8a8a8a", linewidth = 0.30
    )
  }
  if (identical(boundary_shape, "circle")) {
    p_space <- p_space + ggplot2::geom_path(
      data = x$niche_boundary,
      ggplot2::aes(x = axis1, y = axis2,
                   linetype = "Fitted radial boundary"),
      colour = "#333333", linewidth = 0.35
    )
  } else if (identical(boundary_shape, "hull") && nrow(x$current_hull) > 2L) {
    p_space <- p_space + ggplot2::geom_polygon(
      data = x$current_hull,
      ggplot2::aes(x = axis1, y = axis2,
                   linetype = "Current reference climate envelope"),
      fill = "#d8d8d8", colour = "#333333", alpha = 0.72,
      linewidth = 0.30
    )
  }
  if (show_hulls) {
    if (!identical(boundary_shape, "hull") && nrow(x$current_hull) > 2L) {
      p_space <- p_space + ggplot2::geom_polygon(
        data = x$current_hull,
        ggplot2::aes(x = axis1, y = axis2,
                     linetype = "Current reference climate envelope"),
        fill = "#d8d8d8", colour = "#333333", alpha = 0.72,
        linewidth = 0.30
      )
    }
    p_space <- p_space + ggplot2::geom_polygon(
      data = x$future_hull,
      ggplot2::aes(x = axis1, y = axis2,
                   linetype = "Future climate projection of analysed cells"),
      fill = NA, colour = "#1b6f6a", linewidth = 0.40
    )
  }
  if (show_center) {
    p_space <- p_space + ggplot2::geom_point(
      data = x$center,
      ggplot2::aes(x = axis1, y = axis2),
      shape = 21, colour = "black", fill = "white", size = 2.4,
      stroke = 0.55
    )
  }
  p_space <- p_space +
    ggplot2::scale_linetype_manual(
      values = c(
        "Current reference climate envelope" = "solid",
        "Fitted radial boundary" = "dashed",
        "Future climate projection of analysed cells" = "solid"
      ),
      breaks = c(
        "Current reference climate envelope",
        "Fitted radial boundary",
        "Future climate projection of analysed cells"
      ),
      name = "Climate-space envelopes"
    ) +
    ggplot2::guides(
      linetype = if (show_boundary_label) {
        ggplot2::guide_legend(order = 2)
      } else {
        "none"
      }
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "(a) Reduced climate-space projection",
      x = "Reduced niche climate axis 1",
      y = "Reduced niche climate axis 2"
    ) +
    .climniche_theme(base_size = 8.2) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major = ggplot2::element_line(colour = "#e1e1e1",
                                               linewidth = 0.15)
    )

  p_radial <- ggplot2::ggplot(
    cells,
    ggplot2::aes(x = climate_change_amount,
                 y = niche_distance_change,
                 colour = class)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25,
                        linetype = 2, colour = "#8a8a8a") +
    ggplot2::geom_point(size = 0.48, alpha = 0.22) +
    ggplot2::scale_colour_manual(
      values = cell_cols,
      labels = .class_labels(),
      drop = TRUE
    ) +
    ggplot2::labs(
      title = "(b) Displacement and distance shift",
      x = "Climatic Displacement",
      y = "Niche Distance Shift",
      colour = "Derived exposure class"
    ) +
    .climniche_theme(base_size = 8.2)

  p_boundary <- ggplot2::ggplot(
    cells,
    ggplot2::aes(x = climate_reconfiguration,
                 y = niche_boundary_exceedance,
                 colour = class)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25,
                        linetype = 2, colour = "#8a8a8a") +
    ggplot2::geom_point(size = 0.48, alpha = 0.22) +
    ggplot2::scale_colour_manual(
      values = cell_cols,
      labels = .class_labels(),
      drop = TRUE
    ) +
    ggplot2::labs(
      title = "(c) Reconfiguration and boundary exceedance",
      x = "Climatic Reconfiguration",
      y = "Niche Boundary Exceedance",
      colour = "Derived exposure class"
    ) +
    .climniche_theme(base_size = 8.2) +
    ggplot2::theme(legend.position = "none")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    out <- p_space | (p_radial / p_boundary)
    out <- out + patchwork::plot_layout(widths = c(1.05, 1.10))
    if (!is.null(title)) {
      out <- out + patchwork::plot_annotation(title = title)
    }
    return(out)
  }
  list(projection = p_space, displacement = p_radial,
       reconfiguration = p_boundary)
}

#' Plot a niche climate exposure diagram
#'
#' @param x A fitted climniche object or data returned by
#'   `climniche_diagram_data()`.
#' @param scope `"current"` for current reference cells; `"all"` for all
#'   evaluated cells.
#' @param type `"summary"` draws mean current-to-future displacements for
#'   exposure classes; `"sample"` draws sampled cell displacements and future
#'   points.
#' @param summary_layout Layout used when `type = "summary"`. `"metrics"`
#'   draws the climate-space projection and two metric planes; `"ordination"`
#'   draws class mean arrows on the climate-space projection.
#' @param max_arrows Maximum number of current to future arrows to draw when
#'   `type = "sample"`.
#' @param seed Random seed used when subsampling arrows.
#' @param show_reference Logical; draw the full analysed climate-space domain.
#' @param show_hulls Logical; draw current and future climate-space envelopes.
#' @param boundary_shape Boundary display. `"hull"` draws the current reference
#'   climate envelope, `"circle"` draws a constant niche-distance boundary, and
#'   `"none"` suppresses the boundary.
#' @param show_boundary_label Logical; add envelope explanations below the
#'   exposure-class legend.
#' @param show_points Logical; draw future points when `type = "sample"`.
#' @param show_startpoints Logical; draw class mean current positions when
#'   `type = "summary"`.
#' @param show_endpoints Logical; draw class mean future positions when
#'   `type = "summary"`.
#' @param show_center Logical; mark the realised niche centre.
#' @param show_variables Logical; draw environmental variable directions.
#' @param variable_labels Optional named vector replacing variable labels.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#' @export
plot_climniche_diagram <- function(x, scope = c("current", "all"),
                                   type = c("summary", "sample"),
                                   summary_layout = c("metrics", "ordination"),
                                   max_arrows = 350L, seed = 1L,
                                   show_reference = FALSE,
                                   show_hulls = TRUE,
                                   boundary_shape = c("hull", "circle", "none"),
                                   show_boundary_label = TRUE,
                                   show_points = NULL,
                                   show_startpoints = FALSE,
                                   show_endpoints = FALSE,
                                   show_center = TRUE,
                                   show_variables = FALSE,
                                   variable_labels = NULL,
                                   title = NULL) {
  .need_ggplot2()
  type <- match.arg(type)
  summary_layout <- match.arg(summary_layout)
  boundary_shape <- match.arg(boundary_shape)
  if (!inherits(x, "climniche_diagram_data")) {
    x <- climniche_diagram_data(x, scope = scope, max_arrows = max_arrows,
                                seed = seed)
  }
  if (is.null(show_points)) {
    show_points <- identical(type, "sample")
  }
  if (!is.null(variable_labels)) {
    idx <- match(x$variables$variable, names(variable_labels))
    replace <- !is.na(idx)
    x$variables$variable[replace] <- unname(variable_labels[idx[replace]])
  }
  if (identical(type, "summary") && identical(summary_layout, "metrics")) {
    return(.plot_climniche_metric_diagram(
      x = x,
      show_reference = show_reference,
      show_hulls = show_hulls,
      boundary_shape = boundary_shape,
      show_boundary_label = show_boundary_label,
      show_center = show_center,
      title = title,
      max_points = max_arrows * 20L,
      seed = seed
    ))
  }
  cell_cols <- .class_colours()
  p <- ggplot2::ggplot()
  if (show_reference && nrow(x$reference_hull) > 2L) {
    p <- p + ggplot2::geom_polygon(
      data = x$reference_hull,
      ggplot2::aes(x = axis1, y = axis2),
      fill = "#f5f5f5", colour = "#8a8a8a", linewidth = 0.30
    )
  }
  if (identical(boundary_shape, "circle")) {
    p <- p + ggplot2::geom_path(
      data = x$niche_boundary,
      ggplot2::aes(x = axis1, y = axis2,
                   linetype = "Fitted radial boundary"),
      colour = "#333333", linewidth = 0.35
    )
  } else if (identical(boundary_shape, "hull") && nrow(x$current_hull) > 2L) {
    p <- p + ggplot2::geom_polygon(
      data = x$current_hull,
      ggplot2::aes(x = axis1, y = axis2,
                   linetype = "Current reference climate envelope"),
      fill = "#9a9a9a", colour = "#333333", alpha = 0.34,
      linewidth = 0.30
    )
  }
  if (show_hulls) {
    if (!identical(boundary_shape, "hull")) {
      p <- p + ggplot2::geom_polygon(
        data = x$current_hull,
        ggplot2::aes(x = axis1, y = axis2,
                     linetype = "Current reference climate envelope"),
        fill = "#8f8f8f", colour = "#333333", alpha = 0.45,
        linewidth = 0.30
      )
    }
    p <- p + ggplot2::geom_polygon(
      data = x$future_hull,
      ggplot2::aes(x = axis1, y = axis2,
                   linetype = "Future climate projection of analysed cells"),
      fill = NA, colour = "#1b6f6a", linewidth = 0.35
    )
  }
  if (identical(type, "sample")) {
    p <- p + ggplot2::geom_segment(
        data = x$arrows,
        ggplot2::aes(x = current_axis1, y = current_axis2,
                     xend = future_axis1, yend = future_axis2),
        colour = "#4a4a4a", linewidth = 0.22, alpha = 0.32,
        arrow = grid::arrow(length = grid::unit(1.0, "mm"), type = "closed")
      )
  } else {
    curve_values <- c(
      "Limited climate niche change" = -0.35,
      "Closer to current niche" = -0.18,
      "Farther from current niche" = 0.18,
      "Outside current niche boundary" = 0.35,
      "Climatic Reconfiguration with limited Niche Distance Shift" = 0.05
    )
    for (class_name in names(curve_values)) {
      row <- x$class_summary[as.character(x$class_summary$class) == class_name,
                             , drop = FALSE]
      if (!nrow(row)) {
        next
      }
      p <- p + ggplot2::geom_curve(
        data = row,
        ggplot2::aes(x = current_axis1, y = current_axis2,
                     xend = future_axis1, yend = future_axis2,
                     colour = class),
        linewidth = 0.76, alpha = 0.96,
        curvature = unname(curve_values[[class_name]]),
        angle = 90,
        ncp = 8,
        arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed")
      )
    }
    if (show_startpoints) {
      p <- p +
        ggplot2::geom_point(
          data = x$class_summary,
          ggplot2::aes(x = current_axis1, y = current_axis2,
                       shape = "Current class mean"),
          colour = "black",
          fill = "white",
          size = 1.8,
          stroke = 0.42
        )
    }
    if (show_endpoints) {
      p <- p +
        ggplot2::geom_point(
          data = x$class_summary,
          ggplot2::aes(x = future_axis1, y = future_axis2,
                       colour = class,
                       shape = "Future class mean"),
          size = 2.0,
          alpha = 0.9
        )
    }
  }
  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = x$cells,
      ggplot2::aes(x = future_axis1, y = future_axis2, colour = class),
      size = 0.65, alpha = 0.35
    )
  }
  if (show_center) {
    p <- p + ggplot2::geom_point(
      data = x$center,
      ggplot2::aes(x = axis1, y = axis2),
      shape = 21, colour = "black", fill = "white", size = 2.4,
      stroke = 0.55
    )
  }
  p <- p +
    ggplot2::scale_colour_manual(values = cell_cols, drop = TRUE) +
    ggplot2::scale_linetype_manual(
      values = c(
        "Current reference climate envelope" = "solid",
        "Fitted radial boundary" = "dashed",
        "Future climate projection of analysed cells" = "solid"
      ),
      breaks = c(
        "Current reference climate envelope",
        "Fitted radial boundary",
        "Future climate projection of analysed cells"
      ),
      name = "Climate-space envelopes"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        override.aes = list(alpha = 1, size = 1.8, linewidth = 0.5)
      ),
      linetype = if (show_boundary_label) ggplot2::guide_legend(order = 2) else "none"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Reduced niche climate axis 1",
      y = "Reduced niche climate axis 2",
      colour = "Class mean current-to-future\ndisplacement",
      title = title
    ) +
    .climniche_theme(base_size = 8.5) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "#e1e1e1",
                                               linewidth = 0.15),
      legend.key.height = grid::unit(4.6, "mm")
    )

  if (show_startpoints || show_endpoints) {
    p <- p +
      ggplot2::scale_shape_manual(
        values = c("Current class mean" = 21,
                   "Future class mean" = 16),
        name = "Mean positions"
      ) +
      ggplot2::guides(shape = ggplot2::guide_legend(order = 3))
  }

  if (show_variables && nrow(x$variables) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = x$variables,
        ggplot2::aes(x = 0, y = 0, xend = axis1, yend = axis2),
        inherit.aes = FALSE,
        colour = "#2f2f2f", linewidth = 0.3,
        arrow = grid::arrow(length = grid::unit(1.35, "mm"), type = "closed")
      ) +
      ggplot2::geom_text(
        data = x$variables,
        ggplot2::aes(x = axis1, y = axis2, label = variable),
        inherit.aes = FALSE,
        size = 2.3, colour = "#2f2f2f", hjust = -0.08, vjust = -0.08
      )
  }
  p
}
