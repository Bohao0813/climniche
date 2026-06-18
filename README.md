<p align="center" style="margin-top: 28px; margin-bottom: 24px;">
  <img src="man/figures/climniche-hex.svg" width="190" alt="climniche logo" />
</p>

# climniche

`climniche` measures projected climate exposure relative to a species' current
reference niche estimated from occurrence records, range maps, or SDM
suitability weights. It separates the amount of local climatic displacement
from the direction of change relative to the fitted niche centre and the degree
to which future conditions exceed the current empirical niche boundary.

[![R-CMD-check](https://github.com/Bohao0813/climniche/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Bohao0813/climniche/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/climniche)](https://CRAN.R-project.org/package=climniche)
[![pkgdown](https://github.com/Bohao0813/climniche/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Bohao0813/climniche/actions/workflows/pkgdown.yaml)

Website: <https://bohao0813.github.io/climniche/>

![Niche climate exposure concept](man/figures/niche-climate-exposure.svg)

## Installation

Install the CRAN release:

```r
install.packages("climniche")
```

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("Bohao0813/climniche")
```

## Decomposition

`climniche` reports four quantities from one niche-relative decomposition.
Field names in fitted R objects use snake_case; figures and reports use the
formal names below.

- Climatic Displacement (`climate_change_amount`): distance between current
  and future conditions in standardised, sensitivity-weighted climatic space.
- Niche Distance Shift (`niche_distance_change`): signed change in distance
  from the current realised climatic niche centre.
- Climatic Reconfiguration (`climate_reconfiguration`): non radial component
  of climatic displacement. It is derived from Climatic Displacement and Niche
  Distance Shift rather than estimated independently.
- Niche Boundary Exceedance (`niche_boundary_exceedance`): positive excess
  beyond the chosen weighted quantile of current reference-cell distances from
  the realised niche centre.

## Inputs and outputs

`climniche` accepts extracted environmental matrices and aligned `raster` or
`terra` layers. The current reference cells can be supplied as occurrences,
a range raster, a binary SDM raster, or a continuous SDM suitability raster.
Binary rasters are used as 0/1 reference weights. Continuous suitability values
are used as continuous reference weights; `occupied_threshold` only removes low
values and does not convert higher values to 1.

The outputs are cell-level tables, maps, weighted summaries, model-derived
variable contributions, report text, and figure data. The four reported
quantities are continuous. Two optional descriptors summarise whether future
conditions move toward or away from the realised niche centre and whether the
chosen niche boundary is exceeded. The older five-class `classification` field
and `plot_climniche_classes()` remain available for existing workflows.

## Basic use

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

climniche_summary(fit)
climniche_report(fit, species = "example species")
plot_climniche_summary_figure(fit)
```

For spatial data, use `fit_climniche_raster()` or `fit_climniche_terra()`.
Both functions accept binary and continuous reference rasters. `domain` limits
the calculation to the study region, whether it is a terrestrial range, an
island, a catchment, a marine region, or a modelled accessible area. Map
functions accept `study_region` as an `sf`, `Spatial`, `SpatVector`, or x-y
boundary data frame, so terrestrial study boundaries can be drawn without
changing the analysis. Thresholds used by the optional descriptors or legacy
classification remain user-settable.

## Worked example

The [Examples](https://bohao0813.github.io/climniche/articles/climniche-examples.html)
page uses a Mediterranean European anchovy case study with Bio-ORACLE v3 marine
climate layers, prepared future projections, and a continuous SDM suitability
map used as reference weights.

## Contributor

Contributions are welcome. Please make a pull request or contact
bohao.he@polimi.it.
