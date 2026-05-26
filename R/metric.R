#' Build a sensitivity-weighted niche metric
#'
#' @param sensitivity Numeric vector of climate-variable sensitivity weights.
#' @param cnfa Optional CENFA `cnfa` object or list with `sf`, `co`, and `eig`.
#' @param type Metric type. `"diag"` uses variable-level sensitivity weights.
#'   `"factor"` uses a factor metric when CNFA factor coordinates are available.
#'
#' @return A positive semi-definite matrix.
#' @export
niche_metric <- function(sensitivity = NULL, cnfa = NULL,
                         type = c("diag", "factor")) {
  type <- match.arg(type)
  if (is.null(sensitivity)) {
    sensitivity <- .extract_slot(cnfa, "sf")
  }
  if (is.null(sensitivity)) {
    stop("Provide sensitivity or a cnfa object with an sf slot.", call. = FALSE)
  }
  sensitivity <- as.numeric(sensitivity)
  if (any(!is.finite(sensitivity)) || any(sensitivity < 0)) {
    stop("sensitivity must contain finite non-negative values.", call. = FALSE)
  }
  if (all(sensitivity == 0)) {
    stop("At least one sensitivity weight must be positive.", call. = FALSE)
  }

  if (type == "factor") {
    U <- .extract_slot(cnfa, "co")
    rho <- .extract_slot(cnfa, "eig")
    if (!is.null(U) && !is.null(rho)) {
      U <- as.matrix(U)
      rho <- as.numeric(rho)
      if (any(!is.finite(U))) {
        stop("CNFA factor coordinates must be finite.", call. = FALSE)
      }
      if (any(!is.finite(rho)) || any(rho < 0) || mean(rho) <= 0) {
        stop("CNFA factor eigenvalues must be finite non-negative values with positive mean.",
             call. = FALSE)
      }
      if (nrow(U) == length(sensitivity) && ncol(U) == length(rho)) {
        rho <- rho / mean(rho)
        A <- U %*% diag(rho, nrow = length(rho)) %*% t(U)
        return(.validate_niche_metric(A, length(sensitivity), "A"))
      }
      warning("CNFA factor dimensions are inconsistent; falling back to diag metric.",
              call. = FALSE)
    } else {
      warning("CNFA factor slots not available; falling back to diag metric.",
              call. = FALSE)
    }
  }

  w <- sensitivity / mean(sensitivity)
  .validate_niche_metric(diag(w, nrow = length(w)), length(w), "A")
}
