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

`climniche` quantifies niche climate exposure: the amount and direction of
projected climate change at cells currently associated with a taxon, measured
relative to its realised climatic niche.

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

## Metrics

`climniche` separates four cell-level metrics. Field names in fitted R objects
use snake_case; figures and reports use the metric names below.

- Climatic Displacement (`climate_change_amount`): total sensitivity weighted
  climatic displacement.
- Niche Distance Shift (`niche_distance_change`): signed change in distance
  from the current realised climatic niche centre.
- Climatic Reconfiguration (`climate_reconfiguration`): non radial component
  of climatic displacement not captured by change in distance to the niche
  centre.
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

The outputs are cell-level tables, maps, exposure classes, variable
contributions, report text, and figure data.

The exposure classes are derived from the four continuous metrics. They are not
one class per metric. Niche Distance Shift is signed, so it separates movement
toward and away from the realised niche centre. A separate low-change class
identifies cells with limited Climatic Displacement, limited Climatic
Reconfiguration and little Niche Distance Shift. Niche Boundary Exceedance is
treated as an overriding class when future conditions exceed the empirical
boundary.

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
Both functions accept binary and continuous reference rasters. `domain` can be
used to restrict the calculation to a study area such as a marine region or a
modelled accessible area. Classification thresholds can be set directly through
arguments such as `tolerance`, `stable_climate_change`,
`stable_reconfiguration`, `boundary_exceedance_tolerance`, and
`conflict_ratio`.

## Worked example

The [Examples](https://bohao0813.github.io/climniche/articles/climniche-examples.html)
page uses a Mediterranean European anchovy case study with Bio-ORACLE v3 marine
climate layers. It shows VIF-screened multivariate predictors, an SDM
suitability threshold selected by maximum test-set TSS, continuous suitability
weights in the climniche fit, maps of the four metrics, and the derived exposure
classes.

## Contributor

We are welcome any helps! Please make a pull request or reach out to
bohao.he@polimi.it if you want to make any contribution.
