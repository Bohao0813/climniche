# climniche

`climniche` assesses niche climate exposure: projected climate change
interpreted relative to the climate conditions a species currently occupies.
The package can use climate niche factor analysis outputs, but it also works
directly with current and future environmental matrices extracted from
terrestrial or marine rasters.

## Core idea

When climate changes within a species' current distribution, ecologists often
need to know more than whether exposure is large. `climniche` asks whether future
climate at current occurrence locations becomes more similar to, less similar
to, or outside the climate conditions currently associated with the species.

The package compares current and future climate in climate space weighted by
species sensitivity and returns four primary outputs:

- climate change amount: how much local climate changes;
- niche convergence or divergence: whether future climate becomes closer to or
  farther from the species' current realised climatic niche;
- composition change: how much climate changes around the current niche centre
  without changing the overall niche distance;
- niche boundary exceedance: whether future climate exceeds an empirical
  boundary of the species' current climate conditions.

The package also reports a mixed variable response flag when climate variables
contribute in opposite directions.

The package returns maps, cell-level classes, summary tables and report text
that can be used directly or combined with demographic, connectivity or
conservation workflows.

## Minimal example

```r
library(climniche)

sim <- simulate_climniche()

fit <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity
)

report <- climniche_report(fit, species = "example species")
print(report)

plot_climniche_diagram(fit)
plot_climniche_showcase(fit)
plot_climniche_class_summary(fit)
plot_climniche_variable_contribution(fit)
```

`climniche_diagram_data()` returns the data behind the 2D niche climate
diagram: current and future coordinates, class mean arrows, the empirical
niche boundary and class labels.

## Raster workflow

For CENFA compatible examples, use `fit_climniche_raster()` with current and
future `Raster*` climate layers and a presence/range `RasterLayer`. The `occupied`
layer can be binary or continuous. Continuous rasters are thresholded: values
greater than `occupied_threshold` define the current distribution used to
estimate the realised niche. The optional `domain` layer limits where exposure
is evaluated, for example to a marine region, accessible area, management unit
or suitable area predicted by an SDM.

```r
fit <- fit_climniche_raster(
  current = clim_current,
  future = clim_future,
  occupied = occupied_raster,
  domain = analysis_domain_raster,
  center = niche_center,
  sensitivity = variable_sensitivity,
  boundary = 0.95
)

plot_climniche_map(fit, metric = "niche_distance_change", occupied = occupied_raster)
plot_climniche_map(fit, metric = "outside_niche_exceedance", occupied = occupied_raster)
plot_climniche_classes(fit, occupied = occupied_raster)
```

For `terra::SpatRaster` workflows, use `fit_climniche_terra()` with the same
arguments. The map functions accept both `raster::RasterLayer` and one-layer
`terra::SpatRaster` outputs. When `occupied_only = TRUE`, supply the same
`occupied_threshold` used for fitting.

For matrix workflows with extracted SDM values, pass `occupied` as a logical
vector or as the continuous suitability vector and set `occupied_threshold`.

## Report outputs

For report based workflows, `climniche_report()` and
`write_climniche_report()` provide a readable summary of settings, class
proportions, metrics and variable contributions. The full cell level table is
available with `climniche_table()`.

## Interpretation boundaries

The package measures exposure relative to current niche conditions. It does not
estimate persistence, abundance change, dispersal limitation, adaptation or
conservation priority by itself.
