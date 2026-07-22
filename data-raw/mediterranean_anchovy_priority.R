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

# The saved terra pointer is replaced with its on-disk map template.
fit$rasters <- list(
  climate_change_amount = rast(
    file.path(run_dir, "anchovy_climate_change_amount.tif")
  )
)

priority <- climniche_priority(
  fit,
  exposure = "niche_distance_change",
  criterion_name = "Current SDM suitability",
  scope = "current",
  positive_only = TRUE
)
priority_diagnostics <- summary(priority)$diagnostics

priority_summary <- data.frame(
  scope = priority$scope,
  exposure = priority$exposure_label,
  criterion = priority$criterion_label,
  ranked_cells = sum(priority$table$included),
  pareto_fronts = nrow(priority$front_sizes),
  first_front_cells = priority$front_sizes$n[1],
  first_front_fraction = priority_diagnostics$first_front_fraction,
  objective_rank_correlation =
    priority_diagnostics$objective_rank_correlation
)
write.csv(
  priority_summary,
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

priority_colours <- c(
  "#f5f3ed", "#c8ddd2", "#7eb3a2", "#397d70", "#153f3b"
)
priority_map <- as.data.frame(
  priority$rasters$relative_priority,
  xy = TRUE,
  na.rm = TRUE
)
names(priority_map)[3] <- "relative_priority"

degree_label <- function(x, positive, negative) {
  value <- format(abs(x), trim = TRUE, scientific = FALSE)
  hemisphere <- ifelse(x < 0, negative, ifelse(x > 0, positive, ""))
  paste0(value, "\u00b0", hemisphere)
}

plane <- plot_climniche_priority(priority, type = "plane") +
  labs(title = "(a) Priority plane") +
  theme(legend.position = "none")

cell_size <- res(priority$rasters$relative_priority)
map <- ggplot() +
  geom_tile(
    data = priority_map,
    aes(x = x, y = y, fill = relative_priority),
    width = cell_size[1],
    height = cell_size[2]
  ) +
  geom_sf(data = world, fill = "#eeeeee", colour = NA) +
  geom_sf(
    data = med_boundary,
    fill = NA,
    colour = "#333333",
    linewidth = 0.22
  ) +
  scale_fill_gradientn(
    colours = priority_colours,
    limits = c(0, 1),
    name = "Rescaled Pareto depth"
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
  labs(title = "(b) Spatial Pareto depth", x = NULL, y = NULL) +
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

figure <- plane + map + plot_layout(widths = c(0.85, 1.35))
figure_file <- file.path(
  "vignettes", "figures", "anchovy-climniche-priority"
)
ggsave(
  paste0(figure_file, ".png"),
  figure,
  width = 183 / 25.4,
  height = 82 / 25.4,
  dpi = 300
)
ggsave(
  paste0(figure_file, ".pdf"),
  figure,
  width = 183 / 25.4,
  height = 82 / 25.4,
  device = cairo_pdf
)
