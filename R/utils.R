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
