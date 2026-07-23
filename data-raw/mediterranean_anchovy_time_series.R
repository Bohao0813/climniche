r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
invisible(lapply(r_files, source))
suppressPackageStartupMessages({
  library(biooracler)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(terra)
})

sf::sf_use_s2(FALSE)

run_dir <- file.path(
  "output", "mediterranean_anchovy", "20260619_zscore_example"
)
case_dir <- file.path("inst", "extdata", "mediterranean_anchovy")
input_root <- file.path("..", "..", "data-raw")
biooracle_dir <- file.path(input_root, "biooracle_v3")
region_file <- file.path(
  input_root, "marine_regions", "mediterranean_iho_mrgid1905.gpkg"
)
dir.create(biooracle_dir, showWarnings = FALSE, recursive = TRUE)

med_boundary <- st_read(region_file, quiet = TRUE)
med_boundary <- st_transform(med_boundary, 4326)
med_boundary <- st_make_valid(med_boundary)
med_boundary <- st_as_sf(st_union(med_boundary))
med_vect <- terra::vect(med_boundary)
med_bbox <- st_bbox(med_boundary)
region_extent <- terra::ext(
  med_bbox[["xmin"]], med_bbox[["xmax"]],
  med_bbox[["ymin"]], med_bbox[["ymax"]]
)

manifest <- read.csv(
  file.path(case_dir, "anchovy_biooracle_layer_manifest.csv")
)
specs <- manifest[manifest$retained_for_climniche, ]
years <- c(2030, 2050, 2070, 2090)

collapse_time_layers <- function(x) {
  layer_name <- sub("_[0-9]+$", "", names(x))
  x_values <- terra::values(x, mat = TRUE)
  template <- terra::rast(
    nrows = terra::nrow(x),
    ncols = terra::ncol(x),
    xmin = terra::xmin(x),
    xmax = terra::xmax(x),
    ymin = terra::ymin(x),
    ymax = terra::ymax(x),
    crs = terra::crs(x)
  )
  out <- lapply(unique(layer_name), function(name) {
    values <- rowMeans(
      x_values[, layer_name == name, drop = FALSE],
      na.rm = TRUE
    )
    values[is.nan(values)] <- NA_real_
    selected <- terra::setValues(template, values)
    names(selected) <- name
    terra::writeRaster(
      selected,
      tempfile(pattern = paste0(name, "_"), fileext = ".tif"),
      overwrite = TRUE
    )
  })
  do.call(c, out)
}

biooracle_stack <- function(specs, dataset_col, variable_col, time = NULL) {
  pieces <- lapply(unique(specs[[dataset_col]]), function(dataset) {
    use <- specs[[dataset_col]] == dataset
    variables <- unique(specs[[variable_col]][use])
    message("Reading Bio-ORACLE layer ", dataset, ": ",
            paste(variables, collapse = ", "))
    constraints <- list(
      latitude = c(med_bbox[["ymin"]], med_bbox[["ymax"]]),
      longitude = c(med_bbox[["xmin"]], med_bbox[["xmax"]])
    )
    if (!is.null(time)) constraints$time <- c(time, time)

    x <- biooracler::download_layers(
      dataset,
      variables = variables,
      constraints = constraints,
      fmt = "raster",
      directory = biooracle_dir,
      verbose = FALSE
    )
    x <- terra::rast(unique(terra::sources(x)))
    x <- collapse_time_layers(x)
    x[[variables]]
  })
  out <- do.call(c, pieces)
  names(out) <- specs$variable[match(names(out), specs[[variable_col]])]
  out <- terra::crop(out, region_extent)
  terra::crs(out) <- "EPSG:4326"
  terra::mask(out, med_vect)
}

current <- biooracle_stack(
  specs,
  dataset_col = "current_dataset",
  variable_col = "current_variable"
)
future <- lapply(years, function(year) {
  x <- biooracle_stack(
    specs,
    dataset_col = "future_dataset",
    variable_col = "future_variable",
    time = paste0(year, "-01-01T00:00:00Z")
  )
  if (!terra::compareGeom(current, x, stopOnError = FALSE)) {
    x <- terra::resample(x, current, method = "bilinear")
  }
  x
})
names(future) <- paste0("ssp245_", years)

suitability <- terra::rast(file.path(run_dir, "anchovy_sdm_suitability.tif"))
if (!terra::compareGeom(current[[1]], suitability, stopOnError = FALSE)) {
  suitability <- terra::resample(suitability, current[[1]], method = "bilinear")
}

sdm_settings <- read.csv(
  file.path(case_dir, "anchovy_presence_background_sdm_settings.csv")
)
sdm_threshold <- as.numeric(
  sdm_settings$value[sdm_settings$setting == "sdm_threshold"]
)
sensitivity_table <- read.csv(
  file.path(case_dir, "anchovy_climniche_sensitivity_weights.csv")
)
sensitivity <- setNames(
  sensitivity_table$sensitivity_weight,
  sensitivity_table$variable
)

domain <- current[[1]]
terra::values(domain) <- ifelse(
  is.finite(terra::values(domain)),
  1,
  NA_real_
)
names(domain) <- "mediterranean_analysis_domain"

series <- fit_climniche_series(
  current = current,
  future = future,
  time = years,
  scenario = "SSP2-4.5",
  occupied = suitability,
  occupied_threshold = sdm_threshold,
  domain = domain,
  sensitivity = sensitivity
)
saveRDS(series, file.path(run_dir, "anchovy_climniche_series.rds"))

range_summary <- climniche_range_summary(
  series,
  scope = "current",
  area_weight = TRUE
)
persistence <- 2L
departure <- climniche_departure(
  series,
  scope = "current",
  persistence = persistence
)
series_report <- climniche_series_report(
  series,
  species = "European anchovy",
  scope = "current",
  area_weight = TRUE,
  persistence = persistence
)

write.csv(
  range_summary,
  file.path(case_dir, "anchovy_climniche_time_range_summary.csv"),
  row.names = FALSE
)
write.csv(
  series_report$departure_summary,
  file.path(case_dir, "anchovy_climniche_time_departure_summary.csv"),
  row.names = FALSE
)
write.csv(
  series_report$change_rate,
  file.path(case_dir, "anchovy_climniche_time_change_rate.csv"),
  row.names = FALSE
)
write_climniche_series_report(
  series_report,
  file.path(run_dir, "anchovy_climniche_series_report.md")
)

black_theme <- theme(
  text = element_text(colour = "black"),
  plot.title = element_text(face = "bold", colour = "black"),
  axis.text = element_text(colour = "black"),
  axis.title = element_text(colour = "black"),
  axis.line = element_line(colour = "black", linewidth = 0.25),
  axis.ticks = element_line(colour = "black", linewidth = 0.25),
  legend.position = "none"
)

extent_plot <- plot_climniche_time(
  series,
  metric = "exposed_fraction",
  area_weight = TRUE,
  show_models = FALSE
) +
  labs(
    title = "(a) Boundary exceedance fraction",
    x = "Projection year",
    y = "Weighted reference fraction\nbeyond radial niche boundary"
  ) +
  scale_x_continuous(breaks = years) +
  black_theme

severity_plot <- plot_climniche_time(
  series,
  metric = "conditional_relative_exceedance",
  area_weight = TRUE,
  show_models = FALSE
) +
  labs(
    title = "(b) Conditional relative exceedance",
    x = "Projection year",
    y = "Conditional mean relative\nNiche Boundary Exceedance"
  ) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  black_theme

range_plot <- plot_climniche_time(
  series,
  metric = "range_wide_relative_exceedance",
  area_weight = TRUE,
  show_models = FALSE
) +
  labs(
    title = "(c) Range mean relative exceedance",
    x = "Projection year",
    y = "Range mean relative\nNiche Boundary Exceedance"
  ) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  black_theme

time_figure <- extent_plot | severity_plot | range_plot
time_file <- file.path(
  "vignettes", "figures", "anchovy-climniche-through-time"
)
ggsave(
  paste0(time_file, ".png"),
  time_figure,
  width = 183 / 25.4,
  height = 68 / 25.4,
  dpi = 300
)
ggsave(
  paste0(time_file, ".pdf"),
  time_figure,
  width = 183 / 25.4,
  height = 68 / 25.4,
  device = cairo_pdf
)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, 4326)
world <- suppressWarnings(st_crop(world, med_bbox))

degree_label <- function(x, positive, negative) {
  value <- format(abs(x), trim = TRUE, scientific = FALSE)
  hemisphere <- ifelse(x < 0, negative, ifelse(x > 0, positive, ""))
  paste0(value, "\u00b0", hemisphere)
}

reference_cells <- departure$cell
first_departure <- as.numeric(departure$first_persistent_departure)
departure_fraction <- departure$departure_time_fraction

first_raster <- suitability
terra::values(first_raster) <- NA_real_
onset_years <- years[seq_len(length(years) - persistence + 1L)]
first_values <- rep(NA_character_, terra::ncell(first_raster))
first_values[reference_cells] <- ifelse(
  is.na(first_departure),
  "No persistent exceedance",
  as.character(first_departure)
)
first_data <- as.data.frame(first_raster, xy = TRUE, na.rm = FALSE)[, 1:2]
first_data$departure <- factor(
  first_values,
  levels = c(
    as.character(onset_years),
    "No persistent exceedance"
  )
)
first_data <- first_data[!is.na(first_data$departure), ]

fraction_raster <- first_raster
fraction_values <- rep(NA_real_, terra::ncell(fraction_raster))
fraction_values[reference_cells] <- departure_fraction
terra::values(fraction_raster) <- fraction_values
fraction_data <- as.data.frame(fraction_raster, xy = TRUE, na.rm = TRUE)
names(fraction_data)[3] <- "fraction"
cell_size <- terra::res(first_raster)

map_base <- list(
  geom_sf(data = world, fill = "#eeeeee", colour = NA),
  geom_sf(
    data = med_boundary,
    fill = NA,
    colour = "black",
    linewidth = 0.22
  ),
  scale_x_continuous(
    breaks = c(-5, 10, 25),
    labels = function(x) degree_label(x, "E", "W")
  ),
  scale_y_continuous(
    breaks = c(30, 35, 40, 45),
    labels = function(x) degree_label(x, "N", "S")
  ),
  coord_sf(
    xlim = c(med_bbox[["xmin"]], med_bbox[["xmax"]]),
    ylim = c(med_bbox[["ymin"]], med_bbox[["ymax"]]),
    crs = st_crs(4326),
    default_crs = st_crs(4326),
    expand = FALSE
  ),
  theme_classic(base_size = 8.5, base_family = "Arial"),
  theme(
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.25, colour = "black"),
    axis.ticks = element_line(linewidth = 0.25, colour = "black"),
    panel.grid.major = element_line(linewidth = 0.12, colour = "#e2e2e2"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_text(colour = "black"),
    legend.text = element_text(colour = "black"),
    legend.box.margin = margin(t = 2, r = 0, b = 0, l = 0)
  )
)

first_map <- ggplot() +
  geom_tile(
    data = first_data,
    aes(x = x, y = y, fill = departure),
    width = cell_size[1],
    height = cell_size[2]
  ) +
  map_base +
  scale_fill_manual(
    values = c(
      "2030" = "#8c2d04",
      "2050" = "#d94801",
      "2070" = "#f16913",
      "No persistent exceedance" = "#dce8e5"
    ),
    drop = FALSE,
    name = "Projection year"
  ) +
  labs(
    title = "(a) Persistent boundary exceedance onset",
    x = NULL,
    y = NULL
  ) +
  guides(fill = guide_legend(
    nrow = 1,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  ))

fraction_map <- ggplot() +
  geom_tile(
    data = fraction_data,
    aes(x = x, y = y, fill = fraction),
    width = cell_size[1],
    height = cell_size[2]
  ) +
  map_base +
  scale_fill_gradientn(
    colours = c("#f5f3ed", "#b9d8cf", "#559b8c", "#13564f"),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(round(100 * x), "%"),
    name = "Time weighted fraction"
  ) +
  labs(
    title = "(b) Time weighted exceedance fraction",
    x = NULL,
    y = NULL
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(54, "mm"),
    barheight = grid::unit(3.2, "mm")
  ))

spatial_figure <- first_map | fraction_map
spatial_file <- file.path(
  "vignettes", "figures", "anchovy-climniche-time-maps"
)
ggsave(
  paste0(spatial_file, ".png"),
  spatial_figure,
  width = 183 / 25.4,
  height = 91 / 25.4,
  dpi = 300
)
ggsave(
  paste0(spatial_file, ".pdf"),
  spatial_figure,
  width = 183 / 25.4,
  height = 91 / 25.4,
  device = cairo_pdf
)
