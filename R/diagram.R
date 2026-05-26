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

.niche_diagram_rotation <- function(coords, idx) {
  ref <- coords[idx, , drop = FALSE]
  if (ncol(ref) == 1L) {
    return(matrix(c(1, 0), nrow = 1))
  }
  S <- crossprod(ref) / max(1, nrow(ref))
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
#' @param scope `"current"` for current occurrence, range or thresholded SDM
#'   cells; `"all"` for all evaluated cells.
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
  rot <- .niche_diagram_rotation(coords$current, x$occupied)
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
    composition_change = tab$composition_change,
    niche_boundary_exceedance = tab$outside_niche_exceedance,
    class = tab$class,
    mixed_variable_response = tab$mixed_variable_response,
    stringsAsFactors = FALSE
  )
  cells$class <- factor(cells$class, levels = levels(x$classification))
  class_summary <- stats::aggregate(
    cells[, c("current_axis1", "current_axis2", "future_axis1",
              "future_axis2", "climate_change_amount",
              "niche_distance_change", "composition_change",
              "niche_boundary_exceedance")],
    by = list(class = cells$class),
    FUN = mean,
    na.rm = TRUE
  )
  class_counts <- as.data.frame(table(cells$class), stringsAsFactors = FALSE)
  names(class_counts) <- c("class", "count")
  class_summary <- merge(class_summary, class_counts, by = "class",
                         all.x = TRUE, sort = FALSE)
  class_summary$proportion <- class_summary$count / nrow(cells)
  class_summary$class <- factor(class_summary$class,
                                levels = levels(x$classification))

  arrows <- cells
  if (nrow(arrows) > max_arrows) {
    set.seed(seed)
    arrows <- arrows[sample(seq_len(nrow(arrows)), max_arrows), , drop = FALSE]
  }

  reference <- data.frame(axis1 = current2[, 1], axis2 = current2[, 2])
  current_ref <- data.frame(axis1 = current2[x$occupied, 1],
                            axis2 = current2[x$occupied, 2])
  current_radius <- sqrt(rowSums(coords$current[x$occupied, , drop = FALSE]^2))
  current_core <- current_ref[current_radius <= x$boundary_radius, , drop = FALSE]
  if (nrow(current_core) >= 3L) {
    current_ref <- current_core
  }
  future_ref <- data.frame(axis1 = future2[idx, 1], axis2 = future2[idx, 2])
  future_radius <- sqrt(rowSums(coords$future[idx, , drop = FALSE]^2))
  future_cut <- stats::quantile(future_radius, probs = 0.95, na.rm = TRUE,
                                names = FALSE, type = 8)
  future_core <- future_ref[future_radius <= future_cut, , drop = FALSE]
  if (nrow(future_core) >= 3L) {
    future_ref <- future_core
  }

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
      n_cells = nrow(cells),
      n_arrows = nrow(arrows)
    )
  )
  class(out) <- "climniche_diagram_data"
  out
}

#' Plot a niche climate exposure diagram
#'
#' @param x A fitted climniche object or data returned by
#'   `climniche_diagram_data()`.
#' @param scope `"current"` for current occurrence, range or thresholded SDM
#'   cells; `"all"` for all evaluated cells.
#' @param type `"summary"` draws class mean arrows; `"sample"` draws sampled
#'   cell arrows and future points.
#' @param max_arrows Maximum number of current to future arrows to draw when
#'   `type = "sample"`.
#' @param seed Random seed used when subsampling arrows.
#' @param show_reference Logical; draw the full analysed environmental domain.
#' @param show_hulls Logical; draw current and future niche hulls.
#' @param boundary_shape Boundary display. `"hull"` draws the empirical
#'   occupied niche polygon, `"circle"` draws a constant niche distance boundary,
#'   and `"none"` suppresses the boundary.
#' @param show_boundary_label Logical; add reference-area explanations below
#'   the exposure-class legend.
#' @param show_points Logical; draw future points when `type = "sample"`.
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
                                   max_arrows = 350L, seed = 1L,
                                   show_reference = FALSE,
                                   show_hulls = TRUE,
                                   boundary_shape = c("hull", "circle", "none"),
                                   show_boundary_label = TRUE,
                                   show_points = NULL,
                                   show_endpoints = FALSE,
                                   show_center = FALSE,
                                   show_variables = FALSE,
                                   variable_labels = NULL,
                                   title = NULL) {
  .need_ggplot2()
  type <- match.arg(type)
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
                   linetype = "95% current niche distance"),
      colour = "#333333", linewidth = 0.35
    )
  } else if (identical(boundary_shape, "hull") && nrow(x$current_hull) > 2L) {
    p <- p + ggplot2::geom_polygon(
      data = x$current_hull,
      ggplot2::aes(x = axis1, y = axis2,
                   linetype = "95% current niche area"),
      fill = "#9a9a9a", colour = "#333333", alpha = 0.34,
      linewidth = 0.30
    )
  }
  if (show_hulls) {
    if (!identical(boundary_shape, "hull")) {
      p <- p + ggplot2::geom_polygon(
        data = x$current_hull,
        ggplot2::aes(x = axis1, y = axis2,
                     linetype = "current niche area"),
        fill = "#8f8f8f", colour = "#333333", alpha = 0.45,
        linewidth = 0.30
      )
    }
    p <- p + ggplot2::geom_polygon(
      data = x$future_hull,
      ggplot2::aes(x = axis1, y = axis2,
                   linetype = "future climate area"),
      fill = NA, colour = "#1b6f6a", linewidth = 0.35
    )
  }
  if (identical(type, "sample")) {
    p <- p + ggplot2::geom_segment(
        data = x$arrows,
        ggplot2::aes(x = current_axis1, y = current_axis2,
                     xend = future_axis1, yend = future_axis2,
                     colour = class),
        linewidth = 0.22, alpha = 0.32,
        arrow = grid::arrow(length = grid::unit(1.0, "mm"), type = "closed")
      )
  } else {
    p <- p +
      ggplot2::geom_segment(
        data = x$class_summary,
        ggplot2::aes(x = current_axis1, y = current_axis2,
                     xend = future_axis1, yend = future_axis2,
                     colour = class),
        linewidth = 0.75, alpha = 0.95,
        arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed")
      )
    if (show_endpoints) {
      p <- p +
        ggplot2::geom_point(
          data = x$class_summary,
          ggplot2::aes(x = future_axis1, y = future_axis2, colour = class),
          size = 1.8,
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
      shape = 3, colour = "black", size = 2.5, stroke = 0.45
    )
  }
  p <- p +
    ggplot2::scale_colour_manual(values = cell_cols, drop = TRUE) +
    ggplot2::scale_linetype_manual(
      values = c(
        "95% current niche area" = "solid",
        "95% current niche distance" = "dashed",
        "current niche area" = "solid",
        "future climate area" = "solid"
      ),
      breaks = c(
        "95% current niche area",
        "95% current niche distance",
        "current niche area",
        "future climate area"
      ),
      name = "Reference areas"
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
      x = "Niche climate axis 1",
      y = "Niche climate axis 2",
      colour = "Exposure class",
      title = title
    ) +
    .climniche_theme(base_size = 8.5) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "#e1e1e1",
                                               linewidth = 0.15),
      legend.key.height = grid::unit(4.6, "mm")
    )

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
