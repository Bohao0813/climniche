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

.occupied_index <- function(occupied, n, threshold = 0) {
  if (is.null(occupied)) {
    return(seq_len(n))
  }
  if (is.logical(occupied)) {
    if (length(occupied) != n) {
      stop("occupied must have length equal to nrow(current).", call. = FALSE)
    }
    return(which(occupied))
  }
  if (is.numeric(occupied)) {
    if (length(occupied) == n) {
      idx <- which(!is.na(occupied) & occupied > threshold)
      if (!length(idx)) {
        stop("No occupied rows found after applying occupied_threshold.",
             call. = FALSE)
      }
      return(idx)
    }
    occupied <- as.integer(occupied)
    if (any(occupied < 1L | occupied > n)) {
      stop("occupied indices must be between 1 and nrow(current).", call. = FALSE)
    }
    return(unique(occupied))
  }
  stop("occupied must be NULL, a logical vector, or integer indices.", call. = FALSE)
}

.quad_form_rows <- function(x, A) {
  rowSums((x %*% A) * x)
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
