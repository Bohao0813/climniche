#' Build a climatic sensitivity metric
#'
#' @param sensitivity Numeric vector of climate variable sensitivity weights,
#'   used by the diagonal metric. For a factor metric, complete names can be
#'   used to order the CNFA loading rows.
#' @param cnfa Optional CENFA `cnfa` object or compatible list. A diagonal
#'   metric can use `sf`; a factor metric requires `co` and `eig`.
#' @param type Metric type. `"diag"` uses variable-level sensitivity weights.
#'   `"factor"` uses a factor metric when CNFA factor coordinates are available.
#'
#' @return A positive semi-definite matrix.
#'
#' @details
#' For sensitivity weights \eqn{s_j}, the diagonal metric is
#' \deqn{A = \mathrm{diag}(s / \bar{s}).}
#' For a CNFA loading matrix \eqn{U} and factor eigenvalues \eqn{\rho}, the
#' factor metric is
#' \deqn{A = U\,\mathrm{diag}(\rho / \bar{\rho})\,U^{\mathsf{T}}.}
#' The factor metric is constructed from `co` and `eig`; a separate
#' `sensitivity` vector is not combined with it. This quadratic metric is a
#' climniche construction from CENFA components, not a metric returned by
#' CENFA itself.
#'
#' @examples
#' niche_metric(c(temperature = 2, salinity = 1, oxygen = 0.5))
#' @export
niche_metric <- function(sensitivity = NULL, cnfa = NULL,
                         type = c("diag", "factor")) {
  type <- match.arg(type)
  if (type == "factor") {
    U <- .extract_slot(cnfa, "co")
    rho <- .extract_slot(cnfa, "eig")
    if (is.null(U) || is.null(rho)) {
      stop("type = 'factor' requires cnfa co and eig components.",
           call. = FALSE)
    }
    U <- as.matrix(U)
    loading_names <- rownames(U)
    if (!is.null(sensitivity) &&
        (!is.numeric(sensitivity) || any(!is.finite(sensitivity)) ||
         any(sensitivity < 0))) {
      stop("sensitivity must contain finite non-negative values.",
           call. = FALSE)
    }
    sensitivity_names <- names(sensitivity)
    if (.names_are_partial(sensitivity_names) ||
        .names_are_partial(loading_names)) {
      stop("Sensitivity and CNFA loading names must be complete or omitted.",
           call. = FALSE)
    }
    if (!is.null(sensitivity) && length(sensitivity) != nrow(U)) {
      stop("sensitivity must contain one value per CNFA loading row.",
           call. = FALSE)
    }
    if (.names_are_complete(sensitivity_names) &&
        .names_are_complete(loading_names)) {
      .check_unique_names(sensitivity_names, "sensitivity names")
      .check_unique_names(loading_names, "CNFA loading names")
      if (!setequal(sensitivity_names, loading_names)) {
        stop("Sensitivity and CNFA loading names must match.", call. = FALSE)
      }
      U <- U[match(sensitivity_names, loading_names), , drop = FALSE]
    }
    rho_names <- names(rho)
    factor_names <- colnames(U)
    if (.names_are_partial(rho_names) || .names_are_partial(factor_names)) {
      stop("CNFA factor names must be complete or omitted.", call. = FALSE)
    }
    if (.names_are_complete(rho_names) &&
        .names_are_complete(factor_names)) {
      .check_unique_names(rho_names, "CNFA eigenvalue names")
      .check_unique_names(factor_names, "CNFA factor names")
      if (!setequal(rho_names, factor_names)) {
        stop("CNFA factor and eigenvalue names must match.", call. = FALSE)
      }
      rho <- rho[match(factor_names, rho_names)]
    }
    rho <- as.numeric(rho)
    if (any(!is.finite(U))) {
      stop("CNFA factor coordinates must be finite.", call. = FALSE)
    }
    if (any(!is.finite(rho)) || any(rho < 0) || mean(rho) <= 0) {
      stop("CNFA factor eigenvalues must be finite non-negative values with positive mean.",
           call. = FALSE)
    }
    if (ncol(U) != length(rho)) {
      stop("CNFA factor coordinates and eigenvalues have incompatible dimensions.",
           call. = FALSE)
    }
    rho <- rho / mean(rho)
    A <- U %*% diag(rho, nrow = length(rho)) %*% t(U)
    A <- .validate_niche_metric(A, nrow(U), "A")
    if (.names_are_complete(rownames(U))) {
      dimnames(A) <- list(rownames(U), rownames(U))
    }
    return(A)
  }

  if (is.null(sensitivity)) {
    sensitivity <- .extract_slot(cnfa, "sf")
  }
  if (is.null(sensitivity)) {
    stop("Provide sensitivity or a cnfa object with an sf slot.", call. = FALSE)
  }
  sensitivity_names <- names(sensitivity)
  if (.names_are_partial(sensitivity_names)) {
    stop("sensitivity names must be complete or omitted.", call. = FALSE)
  }
  if (.names_are_complete(sensitivity_names)) {
    .check_unique_names(sensitivity_names, "sensitivity names")
  }
  sensitivity <- as.numeric(sensitivity)
  if (any(!is.finite(sensitivity)) || any(sensitivity < 0)) {
    stop("sensitivity must contain finite non-negative values.", call. = FALSE)
  }
  if (all(sensitivity == 0)) {
    stop("At least one sensitivity weight must be positive.", call. = FALSE)
  }
  w <- sensitivity / mean(sensitivity)
  A <- .validate_niche_metric(diag(w, nrow = length(w)), length(w), "A")
  if (.names_are_complete(sensitivity_names)) {
    dimnames(A) <- list(sensitivity_names, sensitivity_names)
  }
  A
}
