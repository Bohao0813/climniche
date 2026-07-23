pkgload::load_all(".")
library(ggplot2)
library(patchwork)
library(sf)
library(terra)

sf::sf_use_s2(FALSE)

run_dir <- file.path(
  "output", "mediterranean_anchovy", "20260619_zscore_example"
)
fit <- readRDS(file.path(run_dir, "anchovy_climniche_fit.rds"))

# Replace the saved terra pointer with the map template on disk.
fit$rasters <- list(
  climate_change_amount = rast(
    file.path(run_dir, "anchovy_climate_change_amount.tif")
  )
)

# Higher SDM suitability with a larger positive niche shift.
positive_shift <- climniche_priority(
  fit,
  exposure = "niche_distance_change",
  criterion_name = "Current SDM suitability",
  scope = "current",
  positive_only = TRUE,
  exposure_direction = "maximize"
)

# Higher SDM suitability with less climatic displacement.
low_displacement <- climniche_priority(
  fit,
  exposure = "climate_change_amount",
  criterion_name = "Current SDM suitability",
  scope = "current",
  positive_only = FALSE,
  exposure_direction = "minimize"
)

priority_summary <- function(x, profile) {
  diagnostics <- summary(x)$diagnostics
  data.frame(
    profile = profile,
    exposure = x$exposure_label,
    exposure_direction = x$exposure_direction,
    criterion = x$criterion_label,
    ranked_cells = sum(x$table$included),
    pareto_fronts = nrow(x$front_sizes),
    first_front_cells = x$front_sizes$n[1],
    first_front_fraction = diagnostics$first_front_fraction,
    objective_rank_correlation = diagnostics$objective_rank_correlation
  )
}

case_summary <- rbind(
  priority_summary(positive_shift, "Positive Niche Distance Shift"),
  priority_summary(
    low_displacement,
    "Low Climatic Displacement"
  )
)
write.csv(
  case_summary,
  file.path(
    "inst", "extdata", "mediterranean_anchovy",
    "anchovy_climniche_priority_summary.csv"
  ),
  row.names = FALSE
)

region_file <- file.path(
  "..", "..", "data-raw", "marine_regions",
  "mediterranean_iho_mrgid1905.gpkg"
)
med_boundary <- st_read(region_file, quiet = TRUE)
med_boundary <- st_transform(med_boundary, 4326)
med_boundary <- st_make_valid(med_boundary)
med_boundary <- st_as_sf(st_union(med_boundary))
med_bbox <- st_bbox(med_boundary)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, 4326)
world <- suppressWarnings(st_crop(world, med_bbox))

positive_shift_colours <- c(
  "#f6f2e8", "#f2ce85", "#e49a48", "#c75b4d", "#7d2431"
)
low_displacement_colours <- c(
  "#f2f5f4", "#c9dfda", "#82b8ae", "#3a8078", "#104a46"
)

degree_label <- function(x, positive, negative) {
  value <- format(abs(x), trim = TRUE, scientific = FALSE)
  hemisphere <- ifelse(x < 0, negative, ifelse(x > 0, positive, ""))
  paste0(value, "\u00b0", hemisphere)
}

priority_map <- function(x, colours, title) {
  raster <- x$rasters$pareto_depth_score
  map_data <- as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  names(map_data)[3] <- "pareto_depth_score"
  cell_size <- res(raster)

  ggplot() +
    geom_tile(
      data = map_data,
      aes(x = x, y = y, fill = pareto_depth_score),
      width = cell_size[1],
      height = cell_size[2]
    ) +
    geom_sf(data = world, fill = "#eeeeee", colour = NA) +
    geom_sf(
      data = med_boundary,
      fill = NA,
      colour = "black",
      linewidth = 0.22
    ) +
    scale_fill_gradientn(
      colours = colours,
      limits = c(0, 1),
      name = "Pareto depth score"
    ) +
    scale_x_continuous(
      breaks = c(-5, 10, 25),
      labels = function(x) degree_label(x, "E", "W")
    ) +
    scale_y_continuous(
      breaks = c(30, 35, 40, 45),
      labels = function(x) degree_label(x, "N", "S")
    ) +
    coord_sf(
      xlim = c(med_bbox[["xmin"]], med_bbox[["xmax"]]),
      ylim = c(med_bbox[["ymin"]], med_bbox[["ymax"]]),
      crs = st_crs(4326),
      default_crs = st_crs(4326),
      expand = FALSE
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 8.5, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      panel.grid.major = element_line(linewidth = 0.12, colour = "#e2e2e2"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(colour = "black"),
      legend.text = element_text(colour = "black"),
      legend.key.width = grid::unit(10, "mm")
    )
}

positive_shift_plane <- plot_climniche_priority(
  positive_shift,
  type = "plane"
) +
  labs(
    title = "(a) Positive Niche Distance Shift: decision plane",
    subtitle = "Higher suitability\nLarger positive Niche Distance Shift"
  ) +
  theme(legend.position = "none")

positive_shift_map <- priority_map(
  positive_shift,
  positive_shift_colours,
  "(b) Positive Niche Distance Shift: Pareto depth"
)

low_displacement_plane <- plot_climniche_priority(
  low_displacement,
  type = "plane"
) +
  labs(
    title = "(c) Low Climatic Displacement: decision plane",
    subtitle = "Higher suitability\nSmaller Climatic Displacement"
  ) +
  theme(legend.position = "none")

low_displacement_map <- priority_map(
  low_displacement,
  low_displacement_colours,
  "(d) Low Climatic Displacement: Pareto depth"
)

row_design <- "
AB
#C
"
positive_shift_row <-
  positive_shift_plane + positive_shift_map + guide_area() +
  plot_layout(
    design = row_design,
    widths = c(0.85, 1.35),
    heights = c(1, 0.13),
    guides = "collect"
  )
low_displacement_row <-
  low_displacement_plane + low_displacement_map + guide_area() +
  plot_layout(
    design = row_design,
    widths = c(0.85, 1.35),
    heights = c(1, 0.13),
    guides = "collect"
  )
figure <- positive_shift_row / low_displacement_row

figure_file <- file.path(
  "vignettes", "figures", "anchovy-climniche-priority"
)
ggsave(
  paste0(figure_file, ".png"),
  figure,
  width = 183 / 25.4,
  height = 150 / 25.4,
  dpi = 300
)
ggsave(
  paste0(figure_file, ".pdf"),
  figure,
  width = 183 / 25.4,
  height = 150 / 25.4,
  device = cairo_pdf
)
