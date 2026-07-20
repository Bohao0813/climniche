.as_numeric_matrix <- function(x, arg = "x") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(arg, " must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(arg, " contains non-finite values; handle NA/Inf before calling climniche.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(arg, " must contain at least one row and one column.",
         call. = FALSE)
  }
  x
}

.names_are_complete <- function(x) {
  !is.null(x) && length(x) > 0L && all(!is.na(x) & nzchar(x))
}

.names_are_partial <- function(x) {
  !is.null(x) && any(!is.na(x) & nzchar(x)) && !.names_are_complete(x)
}

.check_unique_names <- function(x, label) {
  if (anyDuplicated(x)) {
    stop(label, " must be unique.", call. = FALSE)
  }
}

.align_climate_pair <- function(current, future) {
  if (!identical(dim(current), dim(future))) {
    stop("current and future must have identical dimensions.", call. = FALSE)
  }

  current_variables <- colnames(current)
  future_variables <- colnames(future)
  if (.names_are_partial(current_variables) ||
      .names_are_partial(future_variables)) {
    stop("Climate variable names must be complete or omitted.", call. = FALSE)
  }
  if (.names_are_complete(current_variables)) {
    .check_unique_names(current_variables, "current variable names")
  }
  if (.names_are_complete(future_variables)) {
    .check_unique_names(future_variables, "future variable names")
  }
  if (.names_are_complete(current_variables) &&
      .names_are_complete(future_variables)) {
    if (!setequal(current_variables, future_variables)) {
      stop("current and future variable names must match.", call. = FALSE)
    }
    future <- future[, match(current_variables, future_variables), drop = FALSE]
  } else if (.names_are_complete(current_variables)) {
    colnames(future) <- current_variables
  } else if (.names_are_complete(future_variables)) {
    colnames(current) <- future_variables
  }

  current_rows <- rownames(current)
  future_rows <- rownames(future)
  if (.names_are_partial(current_rows) || .names_are_partial(future_rows)) {
    stop("Climate row names must be complete or omitted.", call. = FALSE)
  }
  if (.names_are_complete(current_rows)) {
    .check_unique_names(current_rows, "current row names")
  }
  if (.names_are_complete(future_rows)) {
    .check_unique_names(future_rows, "future row names")
  }
  if (.names_are_complete(current_rows) && .names_are_complete(future_rows)) {
    if (!setequal(current_rows, future_rows)) {
      stop("current and future row names must match.", call. = FALSE)
    }
    future <- future[match(current_rows, future_rows), , drop = FALSE]
  } else if (.names_are_complete(current_rows)) {
    rownames(future) <- current_rows
  } else if (.names_are_complete(future_rows)) {
    rownames(current) <- future_rows
  }

  list(current = current, future = future)
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.check_finite_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(name, " must be a finite number.", call. = FALSE)
  }
  as.numeric(x)
}

.check_open_probability <- function(x, name) {
  x <- .check_finite_scalar(x, name)
  if (x <= 0 || x >= 1) {
    stop(name, " must be between 0 and 1.", call. = FALSE)
  }
  x
}

.check_positive_integer <- function(x, name) {
  x <- .check_finite_scalar(x, name)
  if (x < 1 || x != floor(x)) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(x)
}

.align_metric_inputs <- function(x, center, A) {
  variables <- colnames(x)
  if (.names_are_partial(variables)) {
    stop("Climate variable names must be complete or omitted.",
         call. = FALSE)
  }
  if (.names_are_complete(variables)) {
    .check_unique_names(variables, "climate variable names")
    keep <- rep(TRUE, ncol(x))
    center <- .subset_fit_vector(
      center, keep, ncol(x), "center", variables
    )
    A <- .subset_fit_matrix(A, keep, ncol(x), "A", variables)
  }
  list(center = as.numeric(center), A = as.matrix(A))
}

.clean_reference_weights <- function(weights, threshold = NULL,
                                     arg = "occupied") {
  weights <- as.numeric(weights)
  weights[!is.finite(weights)] <- 0
  if (any(weights < 0)) {
    stop(arg, " contains negative reference weights.", call. = FALSE)
  }
  if (!is.null(threshold)) {
    threshold <- .check_finite_scalar(threshold, "occupied_threshold")
    weights[weights <= threshold] <- 0
  }
  weights
}

.align_reference_values <- function(x, row_names, arg = "occupied") {
  x_names <- names(x)
  if (.names_are_partial(x_names)) {
    stop(arg, " names must be complete or omitted.", call. = FALSE)
  }
  if (!.names_are_complete(x_names)) {
    return(x)
  }
  .check_unique_names(x_names, paste(arg, "names"))
  if (!.names_are_complete(row_names)) {
    return(x)
  }
  if (!setequal(x_names, row_names)) {
    stop(arg, " names must match current row names.", call. = FALSE)
  }
  x[match(row_names, x_names)]
}

.reference_weights <- function(occupied, n, threshold = NULL,
                               row_names = NULL) {
  if (is.null(occupied)) {
    return(rep(1, n))
  }
  if (is.logical(occupied)) {
    if (length(occupied) != n) {
      stop("occupied must have length equal to nrow(current).", call. = FALSE)
    }
    occupied <- .align_reference_values(occupied, row_names)
    weights <- as.numeric(occupied)
    weights[is.na(weights)] <- 0
    return(weights)
  }
  if (is.numeric(occupied)) {
    if (length(occupied) == n) {
      occupied <- .align_reference_values(occupied, row_names)
      return(.clean_reference_weights(occupied, threshold = threshold))
    }
    if (any(!is.finite(occupied)) ||
        any(occupied != floor(occupied)) ||
        any(occupied < 1L | occupied > n)) {
      stop("occupied indices must be between 1 and nrow(current).", call. = FALSE)
    }
    weights <- rep(0, n)
    weights[unique(as.integer(occupied))] <- 1
    return(weights)
  }
  stop("occupied must be NULL, logical, numeric weights, or row indices.",
       call. = FALSE)
}

.positive_reference_indices <- function(weights) {
  idx <- which(weights > 0)
  if (!length(idx)) {
    stop("No reference rows found after applying occupied_threshold.",
         call. = FALSE)
  }
  idx
}

.quad_form_rows <- function(x, A) {
  rowSums((x %*% A) * x)
}

.weighted_mean <- function(x, weights) {
  x <- .as_numeric_matrix(x, "x")
  weights <- .clean_reference_weights(weights, threshold = NULL,
                                      arg = "weights")
  if (length(weights) != nrow(x)) {
    stop("weights must have length equal to nrow(x).", call. = FALSE)
  }
  total <- sum(weights)
  if (!is.finite(total) || total <= 0) {
    stop("weights must contain at least one positive value.", call. = FALSE)
  }
  colSums(x * weights) / total
}

.weighted_mean_vector <- function(x, weights) {
  weights <- .clean_reference_weights(weights, threshold = NULL,
                                      arg = "weights")
  if (length(x) != length(weights)) {
    stop("x and weights must have the same length.", call. = FALSE)
  }
  ok <- is.finite(x) & weights > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(x[ok] * weights[ok]) / sum(weights[ok])
}

.weighted_col_means <- function(x, weights) {
  x <- .as_numeric_matrix(x, "x")
  weights <- .clean_reference_weights(weights, threshold = NULL,
                                      arg = "weights")
  if (length(weights) != nrow(x)) {
    stop("weights must have length equal to nrow(x).", call. = FALSE)
  }
  ok <- weights > 0
  if (!any(ok)) {
    return(rep(NA_real_, ncol(x)))
  }
  colSums(x[ok, , drop = FALSE] * weights[ok]) / sum(weights[ok])
}

.weighted_quantile <- function(x, weights, probs, names = FALSE) {
  weights <- .clean_reference_weights(weights, threshold = NULL,
                                      arg = "weights")
  if (length(x) != length(weights)) {
    stop("x and weights must have the same length.", call. = FALSE)
  }
  probs <- as.numeric(probs)
  if (any(!is.finite(probs)) || any(probs < 0 | probs > 1)) {
    stop("probs must be finite values between 0 and 1.", call. = FALSE)
  }
  ok <- is.finite(x) & weights > 0
  if (!any(ok)) {
    out <- rep(NA_real_, length(probs))
    if (names) names(out) <- paste0(100 * probs, "%")
    return(out)
  }
  x <- x[ok]
  weights <- weights[ok]
  ord <- order(x)
  x <- x[ord]
  weights <- weights[ord]
  cum <- cumsum(weights) / sum(weights)
  out <- vapply(probs, function(prob) {
    index <- which(cum >= prob)[1L]
    if (is.na(index)) index <- length(x)
    x[index]
  }, numeric(1))
  if (names) names(out) <- paste0(100 * probs, "%")
  out
}

.weighted_ecdf_values <- function(reference, values, weights) {
  weights <- .clean_reference_weights(weights, threshold = NULL,
                                      arg = "weights")
  if (length(reference) != length(weights)) {
    stop("reference and weights must have the same length.", call. = FALSE)
  }
  ok <- is.finite(reference) & weights > 0
  if (!any(ok)) {
    return(rep(NA_real_, length(values)))
  }
  ref <- reference[ok]
  w <- weights[ok]
  ord <- order(ref)
  ref <- ref[ord]
  w <- w[ord]
  cum <- cumsum(w) / sum(w)
  idx <- findInterval(values, ref, rightmost.closed = TRUE)
  out <- ifelse(idx > 0, cum[idx], 0)
  out[!is.finite(values)] <- NA_real_
  out
}

.check_probability <- function(x, name) {
  x <- .check_finite_scalar(x, name)
  if (x < 0 || x > 1) {
    stop(name, " must be a finite value between 0 and 1.", call. = FALSE)
  }
  x
}

.fit_reference_weights <- function(x) {
  if (!is.null(x$occupied_weight) &&
      length(x$occupied_weight) == nrow(x$current)) {
    return(x$occupied_weight)
  }
  weights <- rep(0, nrow(x$current))
  weights[x$occupied] <- 1
  weights
}

.fit_exposure_descriptors <- function(x) {
  if (!is.null(x$radial_direction) && !is.null(x$boundary_status)) {
    return(list(
      radial_direction = factor(
        as.character(x$radial_direction),
        levels = .radial_direction_levels()
      ),
      boundary_status = factor(
        as.character(x$boundary_status),
        levels = .boundary_status_levels()
      )
    ))
  }
  settings <- x$descriptor_settings %||% x$threshold_settings
  tolerance <- if (is.null(settings$tolerance)) 0 else settings$tolerance
  boundary_tolerance <- if (is.null(settings$boundary_exceedance_tolerance)) {
    0
  } else {
    settings$boundary_exceedance_tolerance
  }
  .exposure_descriptors(
    niche_distance_change = x$niche_distance_change,
    niche_boundary_exceedance = .fit_metric(x, "niche_boundary_exceedance"),
    tolerance = tolerance,
    tolerance_quantile = settings$tolerance_quantile %||% 0.10,
    boundary_exceedance_tolerance = boundary_tolerance
  )
}

.metric_key <- function(metric) {
  aliases <- c(
    composition_change = "climate_reconfiguration",
    outside_niche_exceedance = "niche_boundary_exceedance"
  )
  if (metric %in% names(aliases)) {
    unname(aliases[[metric]])
  } else {
    metric
  }
}

.fit_metric <- function(x, metric) {
  key <- .metric_key(metric)
  if (!is.null(x[[key]])) {
    return(x[[key]])
  }
  legacy <- c(
    climate_reconfiguration = "composition_change",
    niche_boundary_exceedance = "outside_niche_exceedance"
  )
  if (key %in% names(legacy) && !is.null(x[[legacy[[key]]]])) {
    return(x[[legacy[[key]]]])
  }
  x[[metric]]
}

.validate_niche_metric <- function(A, p, arg = "A",
                                   tol = sqrt(.Machine$double.eps)) {
  A <- as.matrix(A)
  if (!is.numeric(A) || any(!is.finite(A))) {
    stop(arg, " must be a finite numeric matrix.", call. = FALSE)
  }
  if (!identical(dim(A), c(p, p))) {
    stop(arg, " must be a square matrix with dimension ncol(current).",
         call. = FALSE)
  }
  sym_tol <- tol * max(1, max(abs(A)))
  if (max(abs(A - t(A))) > sym_tol) {
    stop(arg, " must be symmetric.", call. = FALSE)
  }
  A <- (A + t(A)) / 2
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  ev_tol <- tol * max(1, max(abs(ev)))
  if (min(ev) < -ev_tol) {
    stop(arg, " must be positive semidefinite.", call. = FALSE)
  }
  if (max(ev) <= ev_tol) {
    stop(arg, " must contain at least one positive dimension.",
         call. = FALSE)
  }
  A
}

.preprocess_climate_pair <- function(current, future, preprocess = TRUE,
                                     correlation = 0.95,
                                     min_sd = 1e-08) {
  p <- ncol(current)
  variable <- colnames(current)
  if (is.null(variable)) {
    variable <- paste0("V", seq_len(p))
    colnames(current) <- variable
    colnames(future) <- variable
  }
  settings <- list(
    enabled = isTRUE(preprocess),
    correlation = correlation,
    min_sd = min_sd
  )
  if (!isTRUE(preprocess)) {
    return(list(
      current = current,
      future = future,
      keep = rep(TRUE, p),
      original_variables = variable,
      retained_variables = variable,
      removed_variables = data.frame(
        variable = character(),
        reason = character(),
        value = numeric(),
        stringsAsFactors = FALSE
      ),
      settings = settings
    ))
  }
  correlation <- .check_finite_scalar(correlation,
                                      "preprocess_correlation")
  if (correlation <= 0 || correlation > 1) {
    stop("preprocess_correlation must be greater than 0 and no larger than 1.",
         call. = FALSE)
  }
  min_sd <- .check_finite_scalar(min_sd, "preprocess_min_sd")
  if (min_sd < 0) {
    stop("preprocess_min_sd must be a finite non-negative number.",
         call. = FALSE)
  }
  settings$correlation <- correlation
  settings$min_sd <- min_sd

  sd_current <- apply(current, 2L, stats::sd)
  keep <- is.finite(sd_current) & sd_current > min_sd
  removed <- data.frame(
    variable = variable[!keep],
    reason = rep("near-zero current variance", sum(!keep)),
    value = sd_current[!keep],
    stringsAsFactors = FALSE
  )
  if (!any(keep)) {
    stop("No climate variables remain after preprocessing.", call. = FALSE)
  }

  active <- which(keep)
  if (length(active) > 1L) {
    cor_current <- stats::cor(current[, active, drop = FALSE])
    cor_current[!is.finite(cor_current)] <- 0
    diag(cor_current) <- 0
    active_keep <- rep(TRUE, length(active))
    repeat {
      ids <- which(active_keep)
      if (length(ids) <= 1L) {
        break
      }
      abs_cor <- abs(cor_current[ids, ids, drop = FALSE])
      diag(abs_cor) <- 0
      max_cor <- max(abs_cor, na.rm = TRUE)
      if (!is.finite(max_cor) || max_cor <= correlation) {
        break
      }
      pairs <- which(upper.tri(abs_cor) & abs_cor == max_cor,
                     arr.ind = TRUE)
      pair_keys <- vapply(seq_len(nrow(pairs)), function(k) {
        paste(sort(variable[active[ids[pairs[k, ]]]]), collapse = "\r")
      }, character(1))
      pair <- pairs[order(pair_keys)[1L], ]
      i <- ids[pair[1]]
      j <- ids[pair[2]]
      others_i <- setdiff(ids, i)
      others_j <- setdiff(ids, j)
      mean_i <- mean(abs(cor_current[i, others_i]), na.rm = TRUE)
      mean_j <- mean(abs(cor_current[j, others_j]), na.rm = TRUE)
      drop_local <- if (is.finite(mean_i) && is.finite(mean_j) &&
                        mean_i != mean_j) {
        if (mean_i > mean_j) i else j
      } else if (variable[active[i]] > variable[active[j]]) {
        i
      } else {
        j
      }
      active_keep[drop_local] <- FALSE
      removed <- rbind(
        removed,
        data.frame(
          variable = variable[active[drop_local]],
          reason = "high current correlation",
          value = max_cor,
          stringsAsFactors = FALSE
        )
      )
    }
    keep[active] <- active_keep
  }
  if (!any(keep)) {
    stop("No climate variables remain after preprocessing.", call. = FALSE)
  }
  list(
    current = current[, keep, drop = FALSE],
    future = future[, keep, drop = FALSE],
    keep = keep,
    original_variables = variable,
    retained_variables = variable[keep],
    removed_variables = removed,
    settings = settings
  )
}

.subset_fit_vector <- function(x, keep, p_old, arg, variables) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x)) {
    stop(arg, " must be numeric.", call. = FALSE)
  }
  retained <- variables[keep]
  x_names <- names(x)
  if (.names_are_partial(x_names)) {
    stop(arg, " names must be complete or omitted.", call. = FALSE)
  }
  if (.names_are_complete(x_names)) {
    .check_unique_names(x_names, paste(arg, "names"))
    target <- if (length(x) == p_old) {
      variables
    } else if (length(x) == sum(keep)) {
      retained
    } else {
      NULL
    }
    if (is.null(target) || !setequal(x_names, target)) {
      stop(arg, " names must match the climate variables.", call. = FALSE)
    }
    x <- x[match(target, x_names)]
  }
  x <- as.numeric(x)
  if (length(x) == p_old) {
    return(x[keep])
  }
  if (length(x) == sum(keep)) {
    return(x)
  }
  stop(arg, " must have length equal to the original or preprocessed number of variables.",
       call. = FALSE)
}

.subset_fit_matrix <- function(x, keep, p_old, arg, variables) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- as.matrix(x)
  p_new <- sum(keep)
  retained <- variables[keep]
  target <- if (identical(dim(x), c(p_old, p_old))) {
    variables
  } else if (identical(dim(x), c(p_new, p_new))) {
    retained
  } else {
    stop(arg, " must match the original or preprocessed number of variables.",
         call. = FALSE)
  }
  row_names <- rownames(x)
  column_names <- colnames(x)
  if (.names_are_partial(row_names) || .names_are_partial(column_names)) {
    stop(arg, " row and column names must be complete or omitted.",
         call. = FALSE)
  }
  if (xor(.names_are_complete(row_names),
          .names_are_complete(column_names))) {
    stop(arg, " must have both row and column names or neither.",
         call. = FALSE)
  }
  if (.names_are_complete(row_names)) {
    .check_unique_names(row_names, paste(arg, "row names"))
    .check_unique_names(column_names, paste(arg, "column names"))
    if (!setequal(row_names, target) || !setequal(column_names, target)) {
      stop(arg, " names must match the climate variables.", call. = FALSE)
    }
    x <- x[match(target, row_names), match(target, column_names), drop = FALSE]
  }
  if (identical(dim(x), c(p_old, p_old))) {
    return(x[keep, keep, drop = FALSE])
  }
  if (identical(dim(x), c(p_new, p_new))) {
    return(x)
  }
  stop(arg, " must match the original or preprocessed number of variables.",
       call. = FALSE)
}

.subset_cnfa_loading <- function(x, keep, p_old, variables) {
  x <- as.matrix(x)
  p_new <- sum(keep)
  retained <- variables[keep]
  target <- if (nrow(x) == p_old) {
    variables
  } else if (nrow(x) == p_new) {
    retained
  } else {
    stop("cnfa co must have rows equal to the original or preprocessed number of variables.",
         call. = FALSE)
  }
  row_names <- rownames(x)
  if (.names_are_partial(row_names)) {
    stop("cnfa co row names must be complete or omitted.", call. = FALSE)
  }
  if (.names_are_complete(row_names)) {
    .check_unique_names(row_names, "cnfa co row names")
    if (!setequal(row_names, target)) {
      stop("cnfa co row names must match the climate variables.",
           call. = FALSE)
    }
    x <- x[match(target, row_names), , drop = FALSE]
  }
  if (nrow(x) == p_old) {
    x <- x[keep, , drop = FALSE]
  }
  x
}

.subset_cnfa_object <- function(cnfa, keep, p_old, variables) {
  if (is.null(cnfa)) {
    return(NULL)
  }
  out <- list(
    mf = .subset_fit_vector(.extract_slot(cnfa, "mf"), keep, p_old,
                            "cnfa mf", variables),
    sf = .subset_fit_vector(.extract_slot(cnfa, "sf"), keep, p_old,
                            "cnfa sf", variables),
    eig = .extract_slot(cnfa, "eig")
  )
  co <- .extract_slot(cnfa, "co")
  if (!is.null(co)) {
    co <- .subset_cnfa_loading(co, keep, p_old, variables)
  }
  out$co <- co
  out
}

.standardize_pair <- function(current, future, scale, global_mean, global_sd) {
  variables <- colnames(current)
  if (!scale) {
    center <- rep(0, ncol(current))
    scale_value <- rep(1, ncol(current))
    names(center) <- variables
    names(scale_value) <- variables
    return(list(current = current, future = future, center = center,
                scale = scale_value))
  }
  if (is.null(global_mean)) {
    global_mean <- colMeans(current)
  }
  if (is.null(global_sd)) {
    global_sd <- apply(current, 2L, stats::sd)
  }
  global_mean <- as.numeric(global_mean)
  global_sd <- as.numeric(global_sd)
  if (length(global_mean) != ncol(current) || any(!is.finite(global_mean))) {
    stop("global_mean must contain one finite value per climate variable.",
         call. = FALSE)
  }
  if (length(global_sd) != ncol(current) || any(!is.finite(global_sd)) ||
      any(global_sd <= 0)) {
    stop("global_sd must contain one finite positive value per climate variable.",
         call. = FALSE)
  }
  names(global_mean) <- variables
  names(global_sd) <- variables
  list(
    current = sweep(sweep(current, 2L, global_mean, "-"), 2L, global_sd, "/"),
    future = sweep(sweep(future, 2L, global_mean, "-"), 2L, global_sd, "/"),
    center = global_mean,
    scale = global_sd
  )
}

.extract_slot <- function(object, name) {
  if (is.null(object)) {
    return(NULL)
  }
  if (methods::is(object, "S4") && name %in% methods::slotNames(object)) {
    return(methods::slot(object, name))
  }
  if (is.list(object) && !is.null(object[[name]])) {
    return(object[[name]])
  }
  NULL
}
