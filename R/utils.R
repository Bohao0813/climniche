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
  x
}

.clean_reference_weights <- function(weights, threshold = NULL,
                                     arg = "occupied") {
  weights <- as.numeric(weights)
  weights[!is.finite(weights)] <- 0
  if (any(weights < 0)) {
    stop(arg, " contains negative reference weights.", call. = FALSE)
  }
  if (!is.null(threshold)) {
    threshold <- as.numeric(threshold)[1]
    if (!is.finite(threshold)) {
      stop("occupied_threshold must be NULL or a finite number.",
           call. = FALSE)
    }
    weights[weights <= threshold] <- 0
  }
  weights
}

.reference_weights <- function(occupied, n, threshold = NULL) {
  if (is.null(occupied)) {
    return(rep(1, n))
  }
  if (is.logical(occupied)) {
    if (length(occupied) != n) {
      stop("occupied must have length equal to nrow(current).", call. = FALSE)
    }
    weights <- as.numeric(occupied)
    weights[is.na(weights)] <- 0
    return(weights)
  }
  if (is.numeric(occupied)) {
    if (length(occupied) == n) {
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
  if (length(unique(weights)) == 1L) {
    return(as.numeric(stats::quantile(
      x, probs = probs, names = names, na.rm = TRUE, type = 8
    )))
  }
  ord <- order(x)
  x <- x[ord]
  weights <- weights[ord]
  cum <- cumsum(weights) / sum(weights)
  out <- stats::approx(
    x = c(0, cum),
    y = c(x[1], x),
    xout = probs,
    method = "constant",
    f = 1,
    ties = "ordered",
    rule = 2
  )$y
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
  x <- as.numeric(x)[1]
  if (!is.finite(x) || x < 0 || x > 1) {
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
  correlation <- as.numeric(correlation)[1]
  if (!is.finite(correlation) || correlation <= 0 || correlation > 1) {
    stop("preprocess_correlation must be greater than 0 and no larger than 1.",
         call. = FALSE)
  }
  min_sd <- as.numeric(min_sd)[1]
  if (!is.finite(min_sd) || min_sd < 0) {
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
      pair <- which(abs_cor == max_cor, arr.ind = TRUE)[1, ]
      i <- ids[pair[1]]
      j <- ids[pair[2]]
      others_i <- setdiff(ids, i)
      others_j <- setdiff(ids, j)
      mean_i <- mean(abs(cor_current[i, others_i]), na.rm = TRUE)
      mean_j <- mean(abs(cor_current[j, others_j]), na.rm = TRUE)
      drop_local <- if (is.finite(mean_i) && is.finite(mean_j) &&
                        mean_i > mean_j) {
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

.subset_fit_vector <- function(x, keep, p_old, arg) {
  if (is.null(x)) {
    return(NULL)
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

.subset_fit_matrix <- function(x, keep, p_old, arg) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- as.matrix(x)
  p_new <- sum(keep)
  if (identical(dim(x), c(p_old, p_old))) {
    return(x[keep, keep, drop = FALSE])
  }
  if (identical(dim(x), c(p_new, p_new))) {
    return(x)
  }
  stop(arg, " must match the original or preprocessed number of variables.",
       call. = FALSE)
}

.subset_cnfa_object <- function(cnfa, keep, p_old) {
  if (is.null(cnfa)) {
    return(NULL)
  }
  out <- list(
    mf = .subset_fit_vector(.extract_slot(cnfa, "mf"), keep, p_old, "cnfa mf"),
    sf = .subset_fit_vector(.extract_slot(cnfa, "sf"), keep, p_old, "cnfa sf"),
    eig = .extract_slot(cnfa, "eig")
  )
  co <- .extract_slot(cnfa, "co")
  if (!is.null(co)) {
    co <- as.matrix(co)
    if (nrow(co) == p_old) {
      co <- co[keep, , drop = FALSE]
    } else if (nrow(co) != sum(keep)) {
      stop("cnfa co must have rows equal to the original or preprocessed number of variables.",
           call. = FALSE)
    }
  }
  out$co <- co
  out
}

.standardize_pair <- function(current, future, scale, global_mean, global_sd) {
  if (!scale) {
    return(list(current = current, future = future, center = rep(0, ncol(current)),
                scale = rep(1, ncol(current))))
  }
  if (is.null(global_mean)) {
    global_mean <- colMeans(current)
  }
  if (is.null(global_sd)) {
    global_sd <- apply(current, 2L, stats::sd)
  }
  if (any(global_sd <= 0)) {
    stop("All climate variables must have positive standard deviation.", call. = FALSE)
  }
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
