#' Extract a tidy climniche table
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current occurrence/range cells or `"all"` for
#'   all evaluated cells.
#'
#' @return A data frame with one row per evaluated cell.
#' @export
climniche_table <- function(x, scope = c("current", "all")) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  idx <- if (scope == "current") x$occupied else seq_len(nrow(x$current))
  mixed <- if (!is.null(x$mixed_variable_response)) {
    x$mixed_variable_response[idx]
  } else {
    rep(FALSE, length(idx))
  }

  data.frame(
    cell = idx,
    current_niche_distance = x$niche_radius_current[idx],
    future_niche_distance = x$niche_radius_future[idx],
    climate_change_amount = x$climate_change_amount[idx],
    niche_distance_change = x$niche_distance_change[idx],
    composition_change = x$composition_change[idx],
    change_alignment = x$change_alignment[idx],
    outside_niche_exceedance = x$outside_niche_exceedance[idx],
    current_niche_percentile = x$niche_percentile$current[idx],
    future_niche_percentile = x$niche_percentile$future[idx],
    percentile_change = x$niche_percentile$delta[idx],
    class = x$classification[idx],
    mixed_variable_response = mixed,
    stringsAsFactors = FALSE
  )
}

#' Summarise climniche results
#'
#' @param x A fitted `climniche_fit` object.
#' @param scope `"current"` for current occurrence/range cells or `"all"` for
#'   all evaluated cells.
#'
#' @return A one-row data frame with key report metrics.
#' @export
climniche_summary <- function(x, scope = c("current", "all")) {
  tab <- climniche_table(x, scope = scope)
  class_prop <- prop.table(table(tab$class))
  get_prop <- function(label) {
    if (label %in% names(class_prop)) unname(class_prop[[label]]) else 0
  }

  data.frame(
    scope = match.arg(scope),
    n = nrow(tab),
    boundary_quantile = x$boundary_quantile,
    boundary_distance = x$boundary_radius,
    mean_climate_change_amount = mean(tab$climate_change_amount, na.rm = TRUE),
    median_climate_change_amount = stats::median(tab$climate_change_amount,
                                                 na.rm = TRUE),
    q90_climate_change_amount = as.numeric(stats::quantile(
      tab$climate_change_amount, 0.90, na.rm = TRUE, names = FALSE
    )),
    mean_niche_distance_change = mean(tab$niche_distance_change, na.rm = TRUE),
    median_niche_distance_change = stats::median(tab$niche_distance_change,
                                                 na.rm = TRUE),
    prop_niche_convergence = get_prop("closer to current niche"),
    prop_niche_divergence = get_prop("farther from current niche"),
    prop_niche_exceedance = get_prop("outside current niche boundary"),
    prop_reconfiguration = get_prop("changed composition, similar distance"),
    prop_stable = get_prop("little climate niche change"),
    prop_mixed_variable_response = mean(tab$mixed_variable_response,
                                        na.rm = TRUE),
    prop_outside_niche = mean(tab$outside_niche_exceedance > 0, na.rm = TRUE),
    mean_outside_niche_exceedance = mean(tab$outside_niche_exceedance,
                                         na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Build a climniche report
#'
#' @param x A fitted `climniche_fit` object.
#' @param species Optional species name used in printed reports.
#' @param scope `"current"` for current occurrence/range cells or `"all"` for
#'   all evaluated cells.
#' @param top_variables Number of variable contributions to show.
#'
#' @return An object of class `climniche_report`.
#' @export
climniche_report <- function(x, species = NULL, scope = c("current", "all"),
                             top_variables = 5) {
  if (!inherits(x, "climniche_fit")) {
    stop("x must be a fitted climniche object.", call. = FALSE)
  }
  scope <- match.arg(scope)
  tab <- climniche_table(x, scope = scope)
  summ <- climniche_summary(x, scope = scope)
  idx <- if (scope == "current") x$occupied else seq_len(nrow(x$current))

  cls <- as.data.frame(prop.table(table(tab$class)), stringsAsFactors = FALSE)
  names(cls) <- c("class", "proportion")
  cls$count <- as.integer(table(tab$class)[as.character(cls$class)])
  cls <- cls[order(cls$proportion, decreasing = TRUE), , drop = FALSE]

  vals <- colMeans(x$variable_contribution[idx, , drop = FALSE], na.rm = TRUE)
  var_tab <- data.frame(
    variable = names(vals),
    mean_contribution = as.numeric(vals),
    abs_mean_contribution = abs(as.numeric(vals)),
    interpretation = ifelse(vals > 0,
                            "less similar to current niche",
                            "more similar to current niche"),
    stringsAsFactors = FALSE
  )
  var_tab <- var_tab[order(var_tab$abs_mean_contribution, decreasing = TRUE),
                     , drop = FALSE]
  var_tab <- utils::head(var_tab, top_variables)

  dominant <- if (nrow(cls)) as.character(cls$class[1]) else NA_character_
  direction <- if (summ$mean_niche_distance_change > 0) {
    "less similar to the current realised climatic niche on average"
  } else if (summ$mean_niche_distance_change < 0) {
    "more similar to the current realised climatic niche on average"
  } else {
    "similar in average distance to the current realised climatic niche"
  }

  interpretation <- c(
    paste0("The analysed cells are projected to become ", direction, "."),
    paste0("The most frequent climniche class is '", dominant, "'."),
    paste0(round(100 * summ$prop_mixed_variable_response, 1),
           "% of analysed cells show mixed variable contributions without ",
           "crossing the niche boundary."),
    paste0(round(100 * summ$prop_outside_niche, 1),
           "% of analysed cells exceed the empirical niche boundary at q = ",
           summ$boundary_quantile, ".")
  )

  out <- list(
    species = species,
    scope = scope,
    settings = data.frame(
      boundary_quantile = x$boundary_quantile,
      boundary_distance = x$boundary_radius,
      n_variables = ncol(x$current),
      n_cells = nrow(tab),
      stringsAsFactors = FALSE
    ),
    summary = summ,
    class_summary = cls,
    top_variables = var_tab,
    interpretation = interpretation,
    table = tab
  )
  class(out) <- "climniche_report"
  out
}

#' @export
print.climniche_report <- function(x, ...) {
  title <- if (is.null(x$species)) {
    "climniche report"
  } else {
    paste0("climniche report: ", x$species)
  }
  cat(title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")
  cat("Scope: ", x$scope, "\n", sep = "")
  cat("Cells: ", x$settings$n_cells, "\n", sep = "")
  cat("Boundary quantile: ", x$settings$boundary_quantile, "\n\n", sep = "")

  cat("Interpretation\n")
  for (line in x$interpretation) {
    cat("- ", line, "\n", sep = "")
  }

  cat("\nSummary\n")
  print(x$summary, row.names = FALSE)

  cat("\nClass proportions\n")
  print(x$class_summary, row.names = FALSE)

  cat("\nTop variable contributions\n")
  print(x$top_variables, row.names = FALSE)
  invisible(x)
}

#' Write a climniche report to Markdown
#'
#' @param report An object returned by [climniche_report()].
#' @param file Output Markdown file.
#'
#' @return Invisibly returns `file`.
#' @export
write_climniche_report <- function(report, file) {
  if (!inherits(report, "climniche_report")) {
    stop("report must be produced by climniche_report().", call. = FALSE)
  }
  title <- if (is.null(report$species)) {
    "Niche climate exposure report"
  } else {
    paste0("Niche climate exposure report: ",
           report$species)
  }

  fmt_row <- function(dat) {
    paste(utils::capture.output(print(dat, row.names = FALSE)), collapse = "\n")
  }

  lines <- c(
    paste0("# ", title),
    "",
    "## Interpretation",
    paste0("- ", report$interpretation),
    "",
    "## Settings",
    "```text",
    fmt_row(report$settings),
    "```",
    "",
    "## Summary",
    "```text",
    fmt_row(report$summary),
    "```",
    "",
    "## Class Proportions",
    "```text",
    fmt_row(report$class_summary),
    "```",
    "",
    "## Top Variable Contributions",
    "```text",
    fmt_row(report$top_variables),
    "```",
    "",
    "## Notes",
    paste(
      "climniche is an exposure interpretation. It does not estimate population",
      "growth, dispersal, adaptation, biotic interactions or conservation",
      "priority by itself."
    )
  )
  writeLines(lines, con = file)
  invisible(file)
}
