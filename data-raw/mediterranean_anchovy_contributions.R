pkgload::load_all(".")

library(ggplot2)
library(patchwork)
library(sf)
library(terra)

sf_use_s2(FALSE)

run_dir <- file.path(
  "output", "mediterranean_anchovy", "20260619_zscore_example"
)
fit <- readRDS(file.path(run_dir, "anchovy_climniche_fit.rds"))

# Restore a live raster template for the saved fit.
fit$rasters <- list(
  climate_change_amount = rast(
    file.path(run_dir, "anchovy_climate_change_amount.tif")
  )
)

contribution <- climniche_dominant_contribution(fit, scope = "current")

case_dir <- file.path("inst", "extdata", "mediterranean_anchovy")
manifest <- read.csv(file.path(case_dir, "anchovy_biooracle_layer_manifest.csv"))
variable_labels <- setNames(manifest$label, manifest$variable)

contribution_summary <- contribution$summary
contribution_summary$label <- unname(
  variable_labels[contribution_summary$variable]
)
contribution_summary <- contribution_summary[, c(
  "variable", "label", "mean_absolute_share",
  "mean_signed_contribution", "positive_contribution_fraction",
  "dominant_weight_fraction", "mean_dominant_share"
)]
write.csv(
  contribution_summary,
  file.path(case_dir, "anchovy_climniche_dominant_contributions.csv"),
  row.names = FALSE
)

region_file <- file.path(
  "..", "..", "data-raw", "marine_regions",
  "mediterranean_iho_mrgid1905.gpkg"
)
mediterranean <- st_read(region_file, quiet = TRUE)
mediterranean <- st_transform(mediterranean, 4326)
mediterranean <- st_make_valid(mediterranean)
mediterranean <- st_as_sf(st_union(mediterranean))
med_bbox <- st_bbox(mediterranean)

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
land <- st_transform(land, 4326)
land <- suppressWarnings(st_crop(land, med_bbox))

dominant_map <- as.data.frame(
  contribution$rasters$dominant_variable,
  xy = TRUE,
  na.rm = TRUE
)
names(dominant_map)[3] <- "variable_id"
dominant_map$variable <- contribution$lookup$variable[
  match(dominant_map$variable_id, contribution$lookup$id)
]
dominant_map$label <- unname(variable_labels[dominant_map$variable])
dominant_map$label[dominant_map$variable == "Tied"] <- "Tied"

share_map <- as.data.frame(
  contribution$rasters$dominant_share,
  xy = TRUE,
  na.rm = TRUE
)
names(share_map)[3] <- "dominant_share"

palette <- c(
  "Mean temperature" = "#0072B2",
  "Temperature range" = "#D55E00",
  "Salinity range" = "#009E73",
  "Mean pH" = "#CC79A7",
  "Mean current speed" = "#56B4E9",
  "Current speed range" = "#6A3D9A",
  "Tied" = "#777777"
)

legend_variables <- c(
  "temperature_mean", "temperature_range", "salinity_range", "ph_mean",
  "sea_water_speed_mean", "sea_water_speed_range"
)
legend_labels <- unname(variable_labels[legend_variables])
dominant_map$label <- factor(
  dominant_map$label,
  levels = c(unname(variable_labels), "Tied")
)

degree_label <- function(x, positive, negative) {
  value <- format(abs(x), trim = TRUE, scientific = FALSE)
  hemisphere <- ifelse(x < 0, negative, ifelse(x > 0, positive, ""))
  paste0(value, "\u00b0", hemisphere)
}

map_theme <- function() {
  theme_classic(base_size = 8.5, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      panel.grid.major = element_line(linewidth = 0.12, colour = "#e2e2e2"),
      legend.position = "bottom",
      legend.title = element_text(colour = "black"),
      legend.text = element_text(colour = "black")
    )
}

map_coordinates <- list(
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
  )
)

cell_size <- res(contribution$rasters$dominant_variable)

variable_plot <- ggplot() +
  geom_tile(
    data = dominant_map,
    aes(x = x, y = y, fill = label),
    width = cell_size[1],
    height = cell_size[2],
    show.legend = TRUE
  ) +
  geom_sf(data = land, fill = "#eeeeee", colour = NA) +
  geom_sf(
    data = mediterranean,
    fill = NA,
    colour = "#333333",
    linewidth = 0.22
  ) +
  scale_fill_manual(
    values = palette,
    breaks = legend_labels,
    drop = FALSE,
    name = "Climate variable"
  ) +
  map_coordinates +
  labs(title = "(a) Dominant climatic contribution", x = NULL, y = NULL) +
  guides(fill = guide_legend(
    nrow = 3,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0
  )) +
  map_theme()

share_plot <- ggplot() +
  geom_tile(
    data = share_map,
    aes(x = x, y = y, fill = dominant_share),
    width = cell_size[1],
    height = cell_size[2]
  ) +
  geom_sf(data = land, fill = "#eeeeee", colour = NA) +
  geom_sf(
    data = mediterranean,
    fill = NA,
    colour = "#333333",
    linewidth = 0.22
  ) +
  scale_fill_gradientn(
    colours = c("#f7f7f4", "#bfd5cf", "#70a89f", "#2d7069", "#123f3c"),
    limits = c(0, 1),
    name = "Dominance share"
  ) +
  map_coordinates +
  labs(title = "(b) Dominance share", x = NULL, y = NULL) +
  guides(fill = guide_colourbar(title.position = "top")) +
  map_theme() +
  theme(legend.key.width = grid::unit(18, "mm"))

figure <- variable_plot + share_plot + plot_layout(widths = c(1.15, 1))
figure_file <- file.path(
  "vignettes", "figures", "anchovy-climniche-contributions"
)
ggsave(
  paste0(figure_file, ".png"),
  figure,
  width = 183 / 25.4,
  height = 88 / 25.4,
  dpi = 300
)
ggsave(
  paste0(figure_file, ".pdf"),
  figure,
  width = 183 / 25.4,
  height = 88 / 25.4,
  device = cairo_pdf
)
