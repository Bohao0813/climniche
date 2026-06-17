#' Simulate a minimal niche relative exposure example
#'
#' @param n Number of climate cells.
#' @param p Number of climate variables.
#' @param seed Random seed.
#' @param rho Pairwise correlation among simulated climate variables.
#' @param prevalence Proportion of background cells treated as true current
#'   occurrence locations under the virtual niche.
#' @param shift Climatic Displacement imposed in the closer to niche and
#'   farther from niche scenarios.
#'
#' @return A list with current, future_toward, future_away, occupied, center,
#'   sensitivity and A.
#' @export
simulate_climniche <- function(n = 2000, p = 2, seed = 1,
                               rho = 0, prevalence = 0.30,
                               shift = 0.40) {
  set.seed(seed)
  if (rho < -1 / (p - 1) || rho >= 1) {
    stop("rho must define a positive-definite equicorrelation matrix.",
         call. = FALSE)
  }
  if (prevalence <= 0 || prevalence >= 1) {
    stop("prevalence must be between 0 and 1.", call. = FALSE)
  }
  if (shift <= 0) {
    stop("shift must be positive.", call. = FALSE)
  }

  center <- rep(0, p)
  sensitivity <- seq(1, 2, length.out = p)
  A <- niche_metric(sensitivity = sensitivity)

  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  L <- chol(Sigma)
  current <- matrix(stats::rnorm(n * p), ncol = p) %*% L
  colnames(current) <- paste0("bio", seq_len(p))
  psi <- niche_potential(current, center = center, A = A)
  occupied <- psi <= stats::quantile(psi, prevalence)

  radial <- sweep(current, 2L, center, "-")
  radial_norm <- sqrt(pmax(.quad_form_rows(radial, A), .Machine$double.eps))
  unit <- radial / radial_norm
  step <- shift * unit
  future_toward <- current - step
  future_away <- current + step

  list(
    current = current,
    future_toward = future_toward,
    future_away = future_away,
    occupied = occupied,
    center = center,
    sensitivity = sensitivity,
    A = A
  )
}
