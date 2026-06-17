<p align="center" style="margin-top: 28px; margin-bottom: 24px;">
  <img src="man/figures/climniche-hex.svg" width="190" alt="climniche logo" />
</p>

# climniche

[![R-CMD-check](https://github.com/Bohao0813/climniche/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Bohao0813/climniche/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/climniche)](https://CRAN.R-project.org/package=climniche)
[![pkgdown](https://github.com/Bohao0813/climniche/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Bohao0813/climniche/actions/workflows/pkgdown.yaml)

Website: <https://bohao0813.github.io/climniche/>

Pages: [Examples](https://bohao0813.github.io/climniche/articles/climniche-examples.html) |
[Reference](https://bohao0813.github.io/climniche/reference/index.html) |
[News](https://bohao0813.github.io/climniche/news/index.html)

![Niche climate exposure concept](man/figures/niche-climate-exposure.svg)

`climniche` decomposes projected climate change relative to a species' current
realised climatic niche. It turns current and future environmental layers,
together with occurrence, range or SDM support, into niche-relative exposure
maps, summaries and reports.

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

`climniche` reports four related quantities from one niche-relative
decomposition. Field names in fitted R objects use snake_case; figures and
reports use the names below.

- Climatic Displacement (`climate_change_amount`): total sensitivity weighted
  climatic displacement.
- Niche Distance Shift (`niche_distance_change`): signed change in distance
  from the current realised climatic niche centre.
- Climatic Reconfiguration (`climate_reconfiguration`): non radial component
  of climatic displacement. It is derived from Climatic Displacement and Niche
  Distance Shift rather than estimated independently.
- Niche Boundary Exceedance (`niche_boundary_exceedance`): positive excess of
  future niche distance beyond the empirical boundary of the current realised
  climatic niche.

## Inputs and outputs

`climniche` accepts extracted environmental matrices and aligned `raster` or
`terra` layers. The current reference cells can be supplied as occurrences,
a range raster, a binary SDM raster, or a continuous SDM suitability raster.
Binary rasters are used as 0/1 reference weights. Continuous suitability values
are used as continuous reference weights; `occupied_threshold` only removes low
values and does not convert higher values to 1.

The outputs are cell-level tables, maps, weighted summaries, variable
contributions, report text, and figure data. The four reported quantities are
continuous. Two optional descriptors summarise Niche Distance Shift as toward,
limited, or away from the realised niche centre and record whether the empirical
niche boundary is exceeded. These descriptors are reported separately because
they answer different questions. The earlier five-class combination remains
available through `classification` and `plot_climniche_classes()` for existing
workflows, but is not required to interpret the decomposition.

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
plot_climniche_showcase(fit)
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
climate layers. It starts from prepared current and future rasters plus a
continuous suitability map, then shows how `climniche` reports the
niche-relative decomposition and spatial summaries.

## Contributor

Contributions are welcome. Please make a pull request or contact
bohao.he@polimi.it.
