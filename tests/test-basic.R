if (requireNamespace("climniche", quietly = TRUE)) {
  library(climniche)
}

sim <- simulate_climniche(n = 500, seed = 10)

toward <- fit_climniche(
  current = sim$current,
  future = sim$future_toward,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity
)

away <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity
)

stopifnot(inherits(toward, "climniche_fit"))
stopifnot(length(toward$climate_change_amount) == nrow(sim$current))
stopifnot(length(toward$composition_change) == nrow(sim$current))
stopifnot(length(toward$change_alignment) == nrow(sim$current))
stopifnot(all(toward$composition_change >= -sqrt(.Machine$double.eps)))
stopifnot(all(abs(toward$change_alignment[!is.na(toward$change_alignment)]) <=
                1 + sqrt(.Machine$double.eps)))
stopifnot(isTRUE(all.equal(
  toward$climate_change_amount^2,
  toward$niche_distance_change^2 + toward$composition_change^2,
  tolerance = 1e-7
)))
stopifnot(mean(toward$niche_distance_change) < mean(away$niche_distance_change))
stopifnot(sum(away$outside_niche_exceedance > 0) >=
            sum(toward$outside_niche_exceedance > 0))

A_full <- matrix(c(1, 0.25, 0.25, 1.4), 2, 2)
fit_full <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  A = A_full
)
stopifnot(isTRUE(all.equal(rowSums(fit_full$variable_contribution),
                           fit_full$psi_future - fit_full$psi_current,
                           tolerance = 1e-7)))
bad_A <- matrix(c(1, 2, 0, 1), 2, 2)
bad_fit <- try(fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  A = bad_A
), silent = TRUE)
stopifnot(inherits(bad_fit, "try-error"))
stopifnot(is.logical(away$mixed_variable_response))
stopifnot(!("mixed variable response" %in% levels(away$classification)))
diagram <- climniche_diagram_data(away, max_arrows = 50)
stopifnot(inherits(diagram, "climniche_diagram_data"))
stopifnot(nrow(diagram$cells) == length(away$occupied))
stopifnot(nrow(diagram$arrows) <= 50)
stopifnot(all(c("current_axis1", "future_axis1", "class") %in%
                names(diagram$cells)))
showcase <- climniche_showcase_data(away, max_points = 80)
stopifnot(inherits(showcase, "climniche_showcase_data"))
stopifnot(nrow(showcase$plane) <= 80)
stopifnot(all(c("x_mid", "y_mid", "cell_weight") %in%
                names(showcase$plane_bins)))
stopifnot(all(c("metric", "value") %in% names(showcase$metrics)))
stopifnot(all(c("proportion", "count") %in% names(showcase$classes)))
stopifnot(all(c("boundary_quantile", "prop_exceeded") %in%
                names(showcase$boundary)))

if (requireNamespace("raster", quietly = TRUE)) {
  r1 <- raster::raster(nrows = 4, ncols = 4, xmn = 0, xmx = 4,
                       ymn = 0, ymx = 4, crs = "+proj=longlat +datum=WGS84")
  raster::values(r1) <- seq_len(raster::ncell(r1))
  r2 <- r1
  raster::values(r2) <- rev(seq_len(raster::ncell(r2)))
  cur <- raster::stack(r1, r2)
  fut <- raster::stack(r1 + 0.5, r2 + 1)
  occ <- r1
  raster::values(occ) <- c(rep(0, 8), rep(0.7, 4), rep(1, 4))
  domain <- r1
  raster::values(domain) <- c(rep(0, 4), rep(1, 12))

  rf <- fit_climniche_raster(cur, fut, occupied = occ, occupied_threshold = 0.5)
  stopifnot(inherits(rf, "climniche_fit"))
  stopifnot(length(rf$occupied) == 8)
  stopifnot(methods::is(rf$rasters$climate_change_amount, "RasterLayer"))
  rf_domain <- fit_climniche_raster(cur, fut, occupied = occ,
                                    occupied_threshold = 0.5,
                                    domain = domain)
  stopifnot(length(rf_domain$climate_change_amount) == 12)
  stopifnot(length(rf_domain$occupied) == 8)
}

if (requireNamespace("terra", quietly = TRUE)) {
  r1 <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
                    ymin = 0, ymax = 4, crs = "EPSG:4326")
  terra::values(r1) <- seq_len(terra::ncell(r1))
  r2 <- r1
  terra::values(r2) <- rev(seq_len(terra::ncell(r2)))
  cur <- c(r1, r2)
  fut <- c(r1 + 0.5, r2 + 1)
  occ <- r1
  terra::values(occ) <- c(rep(0, 8), rep(0.7, 4), rep(1, 4))
  domain <- r1
  terra::values(domain) <- c(rep(0, 4), rep(1, 12))

  tf <- fit_climniche_terra(cur, fut, occupied = occ, occupied_threshold = 0.5)
  stopifnot(inherits(tf, "climniche_fit"))
  stopifnot(length(tf$occupied) == 8)
  stopifnot(methods::is(tf$rasters$climate_change_amount, "SpatRaster"))
  tf_domain <- fit_climniche_terra(cur, fut, occupied = occ,
                                   occupied_threshold = 0.5,
                                   domain = domain)
  stopifnot(length(tf_domain$climate_change_amount) == 12)
  stopifnot(length(tf_domain$occupied) == 8)
}
