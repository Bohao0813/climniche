r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
invisible(lapply(r_files, source))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(biooracler)
  library(usdm)
  library(maxnet)
})

sf::sf_use_s2(FALSE)

species_name <- "Engraulis encrasicolus"
region_name <- "Mediterranean Sea"

input_root <- file.path("..", "..", "data-raw")
data_dir <- file.path(input_root, "mediterranean_anchovy")
biooracle_dir <- file.path(input_root, "biooracle_v3")
region_file <- file.path(
  input_root,
  "marine_regions",
  "mediterranean_iho_mrgid1905.gpkg"
)
run_id <- Sys.getenv("CLIMNICHE_RUN_ID", format(Sys.time(), "%Y%m%d_%H%M%S"))
out_dir <- file.path("output", "mediterranean_anchovy", run_id)

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(biooracle_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(region_file)) {
  stop("Missing input file: ", region_file, call. = FALSE)
}

med_boundary <- sf::st_read(region_file, quiet = TRUE)
med_boundary <- sf::st_transform(med_boundary, 4326)
med_boundary <- sf::st_make_valid(med_boundary)
med_boundary <- sf::st_as_sf(sf::st_union(med_boundary))
med_bbox <- sf::st_bbox(med_boundary)
region_extent <- terra::ext(med_bbox[["xmin"]], med_bbox[["xmax"]],
                            med_bbox[["ymin"]], med_bbox[["ymax"]])
med_vect <- terra::vect(med_boundary)
region_wkt <- paste0(
  "POLYGON((",
  med_bbox[["xmin"]], " ", med_bbox[["ymin"]], ", ",
  med_bbox[["xmax"]], " ", med_bbox[["ymin"]], ", ",
  med_bbox[["xmax"]], " ", med_bbox[["ymax"]], ", ",
  med_bbox[["xmin"]], " ", med_bbox[["ymax"]], ", ",
  med_bbox[["xmin"]], " ", med_bbox[["ymin"]],
  "))"
)

download_obis_occurrences <- function(file, max_records = 7000) {
  if (file.exists(file)) {
    return(utils::read.csv(file, stringsAsFactors = FALSE))
  }
  base <- "https://api.obis.org/v3/occurrence"
  query <- paste(
    paste0("scientificname=", utils::URLencode(species_name, reserved = TRUE)),
    paste0("geometry=", utils::URLencode(region_wkt, reserved = TRUE)),
    paste0("size=", max_records),
    "fields=decimalLongitude,decimalLatitude,scientificName,eventDate",
    "dropped=false",
    "absence=false",
    sep = "&"
  )
  url <- paste0(base, "?", query)
  message("Downloading OBIS records for ", species_name)
  res <- jsonlite::fromJSON(url, flatten = TRUE)
  occ <- as.data.frame(res$results)
  utils::write.csv(occ, file, row.names = FALSE)
  occ
}

occ_file <- file.path(data_dir, "obis_anchovy_mediterranean.csv")
occ <- download_obis_occurrences(occ_file)
names(occ) <- sub("^decimalLongitude$", "longitude", names(occ))
names(occ) <- sub("^decimalLatitude$", "latitude", names(occ))
occ <- occ[is.finite(occ$longitude) & is.finite(occ$latitude), ,
           drop = FALSE]
occ <- occ[occ$longitude >= med_bbox[["xmin"]] &
             occ$longitude <= med_bbox[["xmax"]] &
             occ$latitude >= med_bbox[["ymin"]] &
             occ$latitude <= med_bbox[["ymax"]], , drop = FALSE]
occ <- unique(occ[, intersect(c("scientificName", "longitude", "latitude",
                                "eventDate"), names(occ))])
occ_sf <- sf::st_as_sf(occ, coords = c("longitude", "latitude"),
                       crs = 4326, remove = FALSE)
inside_med <- lengths(sf::st_intersects(occ_sf, med_boundary)) > 0
occ <- occ[inside_med, , drop = FALSE]
utils::write.csv(occ, file.path(out_dir, "anchovy_obis_presence_records.csv"),
                 row.names = FALSE)

future_scenario <- Sys.getenv("BIOORACLE_SCENARIO", "SSP245")
future_year <- as.integer(Sys.getenv("BIOORACLE_YEAR", "2050"))
future_dataset_tag <- tolower(future_scenario)
climniche_scale <- TRUE
climniche_preprocess <- TRUE
climniche_preprocess_correlation <- 0.95
climniche_preprocess_min_sd <- 1e-08
candidate_specs <- data.frame(
  variable = c(
    "temperature_mean", "temperature_range",
    "salinity_mean", "salinity_range",
    "oxygen_mean", "ph_mean", "chlorophyll_mean",
    "sea_water_speed_mean", "sea_water_speed_range"
  ),
  current_dataset = c(
    "thetao_baseline_2000_2019_depthsurf",
    "thetao_baseline_2000_2019_depthsurf",
    "so_baseline_2000_2019_depthsurf",
    "so_baseline_2000_2019_depthsurf",
    "o2_baseline_2000_2018_depthsurf",
    "ph_baseline_2000_2018_depthsurf",
    "chl_baseline_2000_2018_depthsurf",
    "sws_baseline_2000_2019_depthsurf",
    "sws_baseline_2000_2019_depthsurf"
  ),
  current_variable = c(
    "thetao_mean", "thetao_range", "so_mean", "so_range",
    "o2_mean", "ph_mean", "chl_mean", "sws_mean", "sws_range"
  ),
  label = c(
    "Mean temperature",
    "Temperature range",
    "Mean salinity",
    "Salinity range",
    "Mean dissolved oxygen",
    "Mean pH",
    "Mean chlorophyll",
    "Mean current speed",
    "Current speed range"
  ),
  stringsAsFactors = FALSE
)
candidate_specs$future_dataset <- sub(
  "_baseline_[0-9]{4}_[0-9]{4}",
  paste0("_", future_dataset_tag, "_2020_2100"),
  candidate_specs$current_dataset
)
candidate_specs$future_variable <- candidate_specs$current_variable

var_names <- candidate_specs$variable
var_labels <- stats::setNames(candidate_specs$label, candidate_specs$variable)
climate_labels <- var_labels

collapse_time_layers <- function(r) {
  layer_groups <- sub("_[0-9]+$", "", names(r))
  unique_groups <- unique(layer_groups)
  out <- lapply(unique_groups, function(nm) {
    idx <- which(layer_groups == nm)
    if (length(idx) == 1) {
      x <- r[[idx]]
    } else {
      x <- mean(r[[idx]], na.rm = TRUE)
    }
    names(x) <- nm
    x
  })
  do.call(c, out)
}

download_biooracle_stack <- function(specs, dataset_col, variable_col,
                                     time = NULL) {
  pieces <- list()
  datasets <- unique(specs[[dataset_col]])
  for (dataset_id in datasets) {
    vars <- unique(specs[[variable_col]][specs[[dataset_col]] == dataset_id])
    message("Reading Bio-ORACLE layer ", dataset_id, ": ",
            paste(vars, collapse = ", "))
    constraints <- list(
      latitude = c(unname(med_bbox[["ymin"]]), unname(med_bbox[["ymax"]])),
      longitude = c(unname(med_bbox[["xmin"]]), unname(med_bbox[["xmax"]]))
    )
    if (!is.null(time)) {
      constraints$time <- c(time, time)
    }
    r <- biooracler::download_layers(
      dataset_id,
      variables = vars,
      constraints = constraints,
      fmt = "raster",
      directory = biooracle_dir,
      verbose = FALSE
    )
    r <- terra::rast(r)
    r <- collapse_time_layers(r)
    pieces[[dataset_id]] <- r[[vars]]
  }
  out <- do.call(c, pieces)
  names(out) <- specs$variable[match(names(out), specs[[variable_col]])]
  out <- terra::crop(out, region_extent)
  terra::crs(out) <- "EPSG:4326"
  terra::mask(out, med_vect)
}

future_time <- paste0(future_year, "-01-01T00:00:00Z")
climate_current <- download_biooracle_stack(
  candidate_specs,
  dataset_col = "current_dataset",
  variable_col = "current_variable"
)
climate_future <- download_biooracle_stack(
  candidate_specs,
  dataset_col = "future_dataset",
  variable_col = "future_variable",
  time = future_time
)

if (!terra::compareGeom(climate_current, climate_future,
                        stopOnError = FALSE)) {
  climate_future <- terra::resample(climate_future, climate_current,
                                    method = "bilinear")
}

vif_values <- function(x) {
  x <- as.data.frame(x)
  out <- rep(NA_real_, ncol(x))
  names(out) <- names(x)
  if (ncol(x) < 2) {
    return(out)
  }
  ok <- stats::complete.cases(x)
  if (sum(ok) <= ncol(x) + 2) {
    return(out)
  }
  vif_tab <- usdm::vif(x[ok, , drop = FALSE])
  out[match(vif_tab$Variables, names(out))] <- vif_tab$VIF
  out
}

screen_predictors <- function(x, labels, role, priority = names(x),
                              cor_cutoff = 0.85, max_vif = 5,
                              min_retained = 6, max_cells = 20000) {
  vals <- terra::values(x)
  ok <- stats::complete.cases(vals)
  idx <- which(ok)
  if (length(idx) > max_cells) {
    set.seed(42)
    idx <- sample(idx, max_cells)
  }
  mat <- vals[idx, , drop = FALSE]
  cor_mat <- suppressWarnings(stats::cor(mat, use = "pairwise.complete.obs"))
  priority <- intersect(priority, names(x))
  priority <- c(priority, setdiff(names(x), priority))
  min_retained <- min(min_retained, length(priority))
  label_vec <- labels[priority]
  label_vec[is.na(label_vec)] <- priority[is.na(label_vec)]
  retained <- character()
  notes <- setNames(rep("", length(priority)), priority)
  for (nm in priority) {
    max_cor <- if (length(retained) == 0) {
      0
    } else {
      max(abs(cor_mat[nm, retained]), na.rm = TRUE)
    }
    if (length(retained) == 0 || is.na(max_cor) || max_cor <= cor_cutoff) {
      retained <- c(retained, nm)
      notes[nm] <- "retained"
    } else if (length(retained) < min_retained) {
      retained <- c(retained, nm)
      notes[nm] <- paste0("retained to keep at least ", min_retained,
                          " predictors")
    } else {
      notes[nm] <- paste0("removed: absolute correlation > ", cor_cutoff)
    }
  }
  repeat {
    retained_vif <- vif_values(mat[, retained, drop = FALSE])
    if (length(retained) <= min_retained ||
        all(!is.finite(retained_vif)) ||
        max(retained_vif, na.rm = TRUE) <= max_vif) {
      break
    }
    remove_var <- names(which.max(retained_vif))
    retained <- setdiff(retained, remove_var)
    notes[remove_var] <- paste0("removed: VIF > ", max_vif)
  }
  retained_vif <- vif_values(mat[, retained, drop = FALSE])
  if (length(retained) < min_retained ||
      (any(is.finite(retained_vif)) &&
       max(retained_vif, na.rm = TRUE) > max_vif)) {
    stop("Predictor screening for ", role,
         " did not retain at least ", min_retained,
         " predictors with VIF <= ", max_vif, ".",
         call. = FALSE)
  }
  max_cor_all <- vapply(priority, function(nm) {
    others <- setdiff(priority, nm)
    if (length(others) == 0) {
      return(NA_real_)
    }
    max(abs(cor_mat[nm, others]), na.rm = TRUE)
  }, numeric(1))
  data.frame(
    role = role,
    variable = priority,
    label = unname(label_vec),
    retained = priority %in% retained,
    max_abs_correlation = unname(max_cor_all),
    vif_after_screening = unname(retained_vif[match(priority, names(retained_vif))]),
    note = unname(notes[priority]),
    stringsAsFactors = FALSE
  )
}

sdm_candidates <- climate_current
sdm_labels <- climate_labels
exposure_candidates <- climate_current
exposure_future_candidates <- climate_future
candidate_priority <- candidate_specs$variable
predictor_cor_cutoff <- 0.85
predictor_vif_limit <- 5
min_retained_predictors <- 6

sdm_screen <- screen_predictors(
  sdm_candidates,
  labels = sdm_labels,
  role = "SDM suitability model",
  priority = candidate_priority,
  cor_cutoff = predictor_cor_cutoff,
  max_vif = predictor_vif_limit,
  min_retained = min_retained_predictors
)
exposure_screen <- screen_predictors(
  exposure_candidates,
  labels = climate_labels,
  role = "climniche exposure calculation",
  priority = candidate_priority,
  cor_cutoff = predictor_cor_cutoff,
  max_vif = predictor_vif_limit,
  min_retained = min_retained_predictors
)
predictor_screen <- rbind(sdm_screen, exposure_screen)
utils::write.csv(predictor_screen,
                 file.path(out_dir, "anchovy_predictor_correlation_vif_screen.csv"),
                 row.names = FALSE)
sdm_retained <- sdm_screen$variable[sdm_screen$retained]
exposure_retained <- exposure_screen$variable[exposure_screen$retained]
sdm_current <- sdm_candidates[[sdm_retained]]
current <- exposure_candidates[[exposure_retained]]
future <- exposure_future_candidates[[exposure_retained]]
var_names <- exposure_retained
var_labels <- climate_labels[exposure_retained]

complete <- stats::complete.cases(terra::values(current)) &
  stats::complete.cases(terra::values(future)) &
  stats::complete.cases(terra::values(sdm_current))
complete_cells <- which(complete)
presence_cells <- unique(terra::cellFromXY(
  current[[1]],
  as.matrix(occ[, c("longitude", "latitude")])
))
presence_cells <- intersect(presence_cells[!is.na(presence_cells)],
                            complete_cells)

if (length(presence_cells) < 30) {
  stop("Too few OBIS occurrence cells overlap Mediterranean Bio-ORACLE layers.",
       call. = FALSE)
}

.binary_auc <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  if (!any(observed == 1) || !any(observed == 0)) {
    return(NA_real_)
  }
  ranks <- rank(predicted, ties.method = "average")
  n_pos <- sum(observed == 1)
  n_neg <- sum(observed == 0)
  (sum(ranks[observed == 1]) - n_pos * (n_pos + 1) / 2) /
    (n_pos * n_neg)
}

make_sdm <- function(current, presence_cells, complete_cells,
                     n_background = NULL, test_fraction = 0.30) {
  extract_cell_values <- function(x, cells) {
    out <- terra::extract(x, cells)
    out[, setdiff(names(out), "ID"), drop = FALSE]
  }
  set.seed(42)
  background <- setdiff(complete_cells, presence_cells)
  if (is.null(n_background)) {
    n_background <- length(presence_cells)
  }
  if (length(background) > n_background) {
    background <- sample(background, n_background)
  }
  presence_cells <- sample(presence_cells)
  background <- sample(background)
  n_pres_test <- max(1, floor(length(presence_cells) * test_fraction))
  n_bg_test <- max(1, floor(length(background) * test_fraction))
  pres_test_cells <- presence_cells[seq_len(n_pres_test)]
  bg_test_cells <- background[seq_len(n_bg_test)]
  pres_train_cells <- setdiff(presence_cells, pres_test_cells)
  bg_train_cells <- setdiff(background, bg_test_cells)

  pres_env <- extract_cell_values(current, pres_train_cells)
  bg_env <- extract_cell_values(current, bg_train_cells)
  dat <- rbind(pres_env, bg_env)
  pa <- c(rep(1, nrow(pres_env)), rep(0, nrow(bg_env)))
  ok <- stats::complete.cases(dat)
  dat <- as.data.frame(dat[ok, , drop = FALSE])
  pa <- pa[ok]
  feature_classes <- "lq"
  mod <- maxnet::maxnet(
    p = pa,
    data = dat,
    f = maxnet::maxnet.formula(pa, dat, classes = feature_classes)
  )
  pred <- rep(NA_real_, terra::ncell(current[[1]]))
  vals <- as.data.frame(terra::values(current)[complete_cells, ,
                                                drop = FALSE])
  pred[complete_cells] <- stats::predict(mod, vals, type = "cloglog")

  test_env <- rbind(
    extract_cell_values(current, pres_test_cells),
    extract_cell_values(current, bg_test_cells)
  )
  test_pa <- c(rep(1, length(pres_test_cells)),
               rep(0, length(bg_test_cells)))
  test_ok <- stats::complete.cases(test_env)
  test_env <- as.data.frame(test_env[test_ok, , drop = FALSE])
  test_pa <- test_pa[test_ok]
  test_pred <- stats::predict(mod, test_env, type = "cloglog")
  auc <- .binary_auc(test_pa, test_pred)
  thresholds <- unique(as.numeric(stats::quantile(
    test_pred,
    probs = seq(0, 1, length.out = 501),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )))
  eval <- lapply(thresholds, function(th) {
    predicted <- test_pred >= th
    tp <- sum(predicted & test_pa == 1)
    fn <- sum(!predicted & test_pa == 1)
    tn <- sum(!predicted & test_pa == 0)
    fp <- sum(predicted & test_pa == 0)
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    data.frame(
      threshold = th,
      sensitivity = sensitivity,
      specificity = specificity,
      tss = sensitivity + specificity - 1
    )
  })
  eval <- do.call(rbind, eval)
  eval <- eval[is.finite(eval$tss), , drop = FALSE]
  best <- eval[which.max(eval$tss), , drop = FALSE]
  suitability <- current[[1]]
  terra::values(suitability) <- pred
  names(suitability) <- "sdm_suitability"
  list(
    model = mod,
    suitability = suitability,
    threshold = best$threshold,
    evaluation = best,
    threshold_method = "maximum test-set TSS",
    auc = auc,
    n_presence = length(presence_cells),
    n_background = length(background),
    n_presence_train = length(pres_train_cells),
    n_background_train = length(bg_train_cells),
    n_presence_test = length(pres_test_cells),
    n_background_test = length(bg_test_cells),
    test_fraction = test_fraction
  )
}

sdm <- make_sdm(sdm_current, presence_cells, complete_cells)
suitability_r <- sdm$suitability
domain_r <- current[[1]]
terra::values(domain_r) <- ifelse(complete, 1, NA_real_)
domain_r <- terra::mask(domain_r, med_vect)
names(domain_r) <- "mediterranean_analysis_domain"

current_mat <- terra::values(current)[complete_cells, , drop = FALSE]
presence_idx <- which(complete_cells %in% presence_cells)
presence_mat <- current_mat[presence_idx, , drop = FALSE]
sensitivity <- apply(current_mat, 2, stats::var, na.rm = TRUE) /
  pmax(apply(presence_mat, 2, stats::var, na.rm = TRUE),
       .Machine$double.eps)
sensitivity <- pmax(0.25, pmin(4, sensitivity))
sensitivity <- sensitivity / mean(sensitivity, na.rm = TRUE)
names(sensitivity) <- var_names

fit <- fit_climniche_terra(
  current = current,
  future = future,
  occupied = suitability_r,
  occupied_threshold = sdm$threshold,
  domain = domain_r,
  domain_threshold = 0,
  sensitivity = sensitivity,
  boundary = 0.95,
  scale = climniche_scale,
  preprocess = climniche_preprocess,
  preprocess_correlation = climniche_preprocess_correlation,
  preprocess_min_sd = climniche_preprocess_min_sd
)
saveRDS(fit, file.path(out_dir, "anchovy_climniche_fit.rds"))

setting_value <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(NA_character_)
  }
  as.character(x[[1]])
}

collapse_setting <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return("none")
  }
  paste(x, collapse = ", ")
}

preprocess_removed <- fit$preprocessing$removed_variables
removed_predictors <- if (is.null(preprocess_removed) ||
                          !nrow(preprocess_removed)) {
  character()
} else {
  preprocess_removed$variable
}

layer_manifest <- data.frame(
  species = species_name,
  region = region_name,
  variable = candidate_specs$variable,
  label = candidate_specs$label,
  depth = "surface (depthsurf)",
  current_dataset = candidate_specs$current_dataset,
  current_variable = candidate_specs$current_variable,
  current_time = "Bio-ORACLE v3 baseline time layers averaged",
  future_dataset = candidate_specs$future_dataset,
  future_variable = candidate_specs$future_variable,
  future_scenario = future_scenario,
  future_time = future_time,
  retained_for_sdm = candidate_specs$variable %in% sdm_retained,
  retained_for_climniche = candidate_specs$variable %in% exposure_retained,
  stringsAsFactors = FALSE
)
utils::write.csv(layer_manifest,
                 file.path(out_dir, "anchovy_biooracle_layer_manifest.csv"),
                 row.names = FALSE)

sdm_settings <- data.frame(
  setting = c(
    "species",
    "region",
    "occurrence_source",
    "occurrence_record_count",
    "presence_cell_count",
    "background_definition",
    "background_cell_count",
    "background_to_presence_ratio",
    "test_fraction",
    "presence_train_count",
    "background_train_count",
    "presence_test_count",
    "background_test_count",
    "sdm_algorithm",
    "sdm_output_used_by_climniche",
    "sdm_threshold_method",
    "sdm_threshold",
    "sdm_test_auc",
    "sdm_test_tss",
    "sdm_test_sensitivity",
    "sdm_test_specificity"
  ),
  value = c(
    species_name,
    region_name,
    "OBIS presence records clipped to the Mediterranean Sea",
    as.character(nrow(occ)),
    as.character(sdm$n_presence),
    paste(
      "random background cells from analysable Mediterranean cells",
      "excluding occurrence cells"
    ),
    as.character(sdm$n_background),
    formatC(sdm$n_background / sdm$n_presence, format = "f", digits = 2),
    as.character(sdm$test_fraction),
    as.character(sdm$n_presence_train),
    as.character(sdm$n_background_train),
    as.character(sdm$n_presence_test),
    as.character(sdm$n_background_test),
    "maxnet presence-background model",
    "continuous suitability used as reference weights",
    sdm$threshold_method,
    formatC(sdm$threshold, format = "f", digits = 4),
    formatC(sdm$auc, format = "f", digits = 3),
    formatC(sdm$evaluation$tss, format = "f", digits = 3),
    formatC(sdm$evaluation$sensitivity, format = "f", digits = 3),
    formatC(sdm$evaluation$specificity, format = "f", digits = 3)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(sdm_settings,
                 file.path(out_dir, "anchovy_presence_background_sdm_settings.csv"),
                 row.names = FALSE)

fit_settings <- data.frame(
  setting = c(
    "species",
    "region",
    "predictor_correlation_cutoff",
    "predictor_vif_limit",
    "vif_package",
    "minimum_retained_predictors",
    "sdm_retained_predictor_count",
    "exposure_retained_predictor_count",
    "sdm_retained_predictors",
    "exposure_retained_predictors",
    "standardisation",
    "preprocess",
    "preprocess_correlation",
    "preprocess_min_sd",
    "preprocess_retained_predictor_count",
    "preprocess_removed_predictor_count",
    "preprocess_retained_predictors",
    "preprocess_removed_predictors",
    "boundary_quantile",
    "tolerance",
    "tolerance_quantile",
    "boundary_exceedance_tolerance"
  ),
  value = c(
    species_name,
    region_name,
    as.character(predictor_cor_cutoff),
    as.character(predictor_vif_limit),
    "usdm::vif",
    as.character(min_retained_predictors),
    as.character(length(sdm_retained)),
    as.character(length(exposure_retained)),
    paste(sdm_retained, collapse = ", "),
    paste(exposure_retained, collapse = ", "),
    if (isTRUE(climniche_scale)) {
      "Z-score using current-layer means and standard deviations"
    } else {
      "not standardised"
    },
    as.character(fit$preprocessing$settings$enabled),
    setting_value(fit$preprocessing$settings$correlation),
    setting_value(fit$preprocessing$settings$min_sd),
    as.character(length(fit$preprocessing$retained_variables)),
    as.character(length(removed_predictors)),
    collapse_setting(fit$preprocessing$retained_variables),
    collapse_setting(removed_predictors),
    setting_value(fit$boundary_quantile),
    setting_value(fit$descriptor_settings$tolerance),
    setting_value(fit$descriptor_settings$tolerance_quantile),
    setting_value(fit$descriptor_settings$boundary_exceedance_tolerance)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(fit_settings,
                 file.path(out_dir, "anchovy_climniche_fit_settings.csv"),
                 row.names = FALSE)

report <- climniche_report(fit, species = species_name, scope = "current")
write_climniche_report(report, file.path(out_dir, "anchovy_climniche_report.md"))
utils::write.csv(climniche_summary(fit, scope = "current"),
                 file.path(out_dir, "anchovy_climniche_metric_summary.csv"),
                 row.names = FALSE)
utils::write.csv(report$descriptor_summary,
                 file.path(out_dir, "anchovy_climniche_descriptor_proportions.csv"),
                 row.names = FALSE)
utils::write.csv(report$top_variables,
                 file.path(out_dir, "anchovy_climniche_variable_contributions.csv"),
                 row.names = FALSE)
utils::write.csv(climniche_table(fit, scope = "current"),
                 file.path(out_dir, "anchovy_climniche_cell_table.csv"),
                 row.names = FALSE)

amount_r <- fit$rasters$climate_change_amount
distance_r <- fit$rasters$niche_distance_change
composition_r <- fit$rasters$climate_reconfiguration
exceed_r <- fit$rasters$niche_boundary_exceedance

terra::writeRaster(suitability_r, file.path(out_dir, "anchovy_sdm_suitability.tif"),
                   overwrite = TRUE)
terra::writeRaster(amount_r, file.path(out_dir, "anchovy_climate_change_amount.tif"),
                   overwrite = TRUE)
terra::writeRaster(distance_r, file.path(out_dir, "anchovy_niche_distance_change.tif"),
                   overwrite = TRUE)
terra::writeRaster(composition_r, file.path(out_dir, "anchovy_climate_reconfiguration.tif"),
                   overwrite = TRUE)
terra::writeRaster(exceed_r, file.path(out_dir, "anchovy_boundary_exceedance.tif"),
                   overwrite = TRUE)

map_df <- function(r, name) {
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d) <- c("x", "y", name)
  d
}

suitability_df <- map_df(suitability_r, "value")
suitability_df <- suitability_df[suitability_df$value > sdm$threshold, ,
                                 drop = FALSE]
amount_df <- map_df(amount_r, "value")
distance_df <- map_df(distance_r, "value")
composition_df <- map_df(composition_r, "value")
exceed_df <- map_df(exceed_r, "value")

suitable_mask <- suitability_df[, c("x", "y"), drop = FALSE]
suitable_key <- paste(suitable_mask$x, suitable_mask$y, sep = ":")
mask_suitable_df <- function(d) {
  key <- paste(d$x, d$y, sep = ":")
  d[key %in% suitable_key, , drop = FALSE]
}
amount_df <- mask_suitable_df(amount_df)
distance_df <- mask_suitable_df(distance_df)
composition_df <- mask_suitable_df(composition_df)
exceed_df <- mask_suitable_df(exceed_df)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- sf::st_transform(world, 4326)
region_crop_box <- sf::st_bbox(
  c(xmin = unname(med_bbox[["xmin"]]),
    ymin = unname(med_bbox[["ymin"]]),
    xmax = unname(med_bbox[["xmax"]]),
    ymax = unname(med_bbox[["ymax"]])),
  crs = sf::st_crs(4326)
)
world_plot <- suppressWarnings(sf::st_crop(world, region_crop_box))

format_degree <- function(x, positive, negative) {
  value <- abs(x)
  value <- ifelse(abs(value - round(value)) < 1e-6,
                  as.character(round(value)),
                  formatC(value, format = "f", digits = 1))
  hemi <- ifelse(x < 0, negative, ifelse(x > 0, positive, ""))
  paste0(value, "\u00b0", hemi)
}

format_lon <- function(x) format_degree(x, "E", "W")
format_lat <- function(x) format_degree(x, "N", "S")

theme_map <- function(base_size = 7.0, show_x = TRUE, show_y = TRUE) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 0.2,
                                colour = "black", hjust = 0,
                                margin = margin(b = 0.8)),
      axis.title = element_blank(),
      axis.text = element_text(size = base_size - 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.22, colour = "black"),
      axis.line = element_line(linewidth = 0.22, colour = "black"),
      panel.grid.major = element_line(linewidth = 0.12, colour = "#e2e2e2"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1.0, colour = "black"),
      legend.key.height = grid::unit(1.8, "mm"),
      legend.key.width = grid::unit(5.0, "mm"),
      legend.margin = margin(1.2, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(1.4, "mm"),
      legend.spacing.y = grid::unit(0, "mm"),
      plot.margin = margin(0.4, 0.8, 0.4, 0.8),
      axis.text.x = if (show_x) element_text(colour = "black") else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.22, colour = "black") else element_blank(),
      axis.text.y = if (show_y) element_text(colour = "black") else element_blank(),
      axis.ticks.y = if (show_y) element_line(linewidth = 0.22, colour = "black") else element_blank()
    )
}

guide_metric <- function(bar_width_mm = 32) {
  guide_colourbar(
    direction = "horizontal",
    barwidth = grid::unit(bar_width_mm, "mm"),
    barheight = grid::unit(1.8, "mm"),
    frame.colour = "black",
    ticks.colour = "black"
  )
}

base_map <- function() {
  ggplot() +
    scale_x_continuous(breaks = c(-5, 10, 25),
                       labels = format_lon) +
    scale_y_continuous(breaks = c(30, 35, 40, 45),
                       labels = format_lat) +
    coord_sf(xlim = c(med_bbox[["xmin"]], med_bbox[["xmax"]]),
             ylim = c(med_bbox[["ymin"]], med_bbox[["ymax"]]),
             crs = sf::st_crs(4326),
             default_crs = sf::st_crs(4326),
             expand = FALSE)
}

geography_layers <- function() {
  list(
    geom_sf(data = world_plot, inherit.aes = FALSE,
            fill = "#eeeeee", colour = NA),
    geom_sf(data = med_boundary, inherit.aes = FALSE, fill = NA,
            colour = "#333333", linewidth = 0.22)
  )
}

cell_width <- terra::res(current)[1]
cell_height <- terra::res(current)[2]

p_sdm <- base_map() +
  geom_tile(data = suitability_df, aes(x = x, y = y, fill = value),
            width = cell_width, height = cell_height, alpha = 0.95) +
  geography_layers() +
  scale_fill_gradientn(
    colours = c("#eef7f2", "#b9dfcf", "#5aa68f", "#1e6e66"),
    name = NULL
  ) +
  guides(fill = guide_metric()) +
  labs(title = "Current suitability weights") +
  theme_map(show_x = FALSE, show_y = TRUE)

p_amount <- base_map() +
  geom_tile(data = amount_df, aes(x = x, y = y, fill = value),
            width = cell_width, height = cell_height, alpha = 0.95) +
  geography_layers() +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#d6e6f2", "#91b9d5", "#27658f"),
    limits = c(0, max(amount_df$value, na.rm = TRUE)),
    name = NULL
  ) +
  guides(fill = guide_metric()) +
  labs(title = "(a) Climatic Displacement") +
  theme_map(show_x = FALSE, show_y = TRUE)

distance_limits <- range(distance_df$value, na.rm = TRUE)
p_distance <- base_map() +
  geom_tile(data = distance_df, aes(x = x, y = y, fill = value),
            width = cell_width, height = cell_height, alpha = 0.95) +
  geography_layers() +
  scale_fill_gradientn(
    colours = c("#fff7f3", "#fdd0c7", "#f1695b", "#a33430"),
    limits = distance_limits, oob = scales::squish,
    name = NULL
  ) +
  guides(fill = guide_metric()) +
  labs(title = "(b) Niche Distance Shift") +
  theme_map(show_x = FALSE, show_y = FALSE)

p_composition <- base_map() +
  geom_tile(data = composition_df, aes(x = x, y = y, fill = value),
            width = cell_width, height = cell_height, alpha = 0.95) +
  geography_layers() +
  scale_fill_gradientn(
    colours = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45"),
    limits = c(0, max(composition_df$value, na.rm = TRUE)),
    name = NULL
  ) +
  guides(fill = guide_metric()) +
  labs(title = "(c) Climatic Reconfiguration") +
  theme_map(show_x = TRUE, show_y = TRUE)

p_exceed <- base_map() +
  geom_tile(data = exceed_df, aes(x = x, y = y, fill = value),
            width = cell_width, height = cell_height, alpha = 0.95) +
  geography_layers() +
  scale_fill_gradientn(
    colours = c("white", "#f3dfb8", "#d9942f", "#8a3f20"),
    limits = c(0, max(exceed_df$value, na.rm = TRUE)),
    name = NULL
  ) +
  guides(fill = guide_metric()) +
  labs(title = "(d) Niche Boundary Exceedance") +
  theme_map(show_x = TRUE, show_y = FALSE)

fig_main <- (p_amount | p_distance) / (p_composition | p_exceed) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))

diagram_labels <- var_labels

fig_report <- plot_climniche_summary_figure(
  fit,
  scope = "current",
  plane_bins = 45,
  variable_labels = diagram_labels,
  title = NULL
)

save_plot <- function(plot, file_base, width_mm, height_mm, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  ggplot2::ggsave(paste0(file_base, ".png"), plot,
                  width = w, height = h, dpi = dpi)
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot,
                  width = w, height = h, device = grDevices::cairo_pdf)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(paste0(file_base, ".svg"), plot,
                    width = w, height = h, device = svglite::svglite)
  }
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(paste0(file_base, ".tiff"), plot,
                    width = w, height = h, device = ragg::agg_tiff,
                    dpi = dpi, compression = "lzw")
  }
}

save_plot(fig_main,
          file.path(out_dir, "figure_anchovy_mediterranean_climniche_maps"),
          width_mm = 183, height_mm = 102)
save_plot(fig_report,
          file.path(out_dir, "figure_anchovy_mediterranean_climniche_summary"),
          width_mm = 183, height_mm = 128)

cat("\nOutput directory:", normalizePath(out_dir, winslash = "/"), "\n")
cat("Clean OBIS records inside Mediterranean:", nrow(occ), "\n")
cat("Presence cells:", length(presence_cells), "\n")
cat("SDM threshold:", round(sdm$threshold, 4), "\n")
print(climniche_summary(fit, scope = "current"))
print(report$descriptor_summary)
