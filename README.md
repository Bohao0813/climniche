<p align="center" style="margin-top: 28px; margin-bottom: 24px;">
  <img src="man/figures/climniche-hex.svg" width="190" alt="climniche logo" />
</p>

# climniche

Equal climatic displacement can move two sites in opposite directions relative
to a species' realised niche. `climniche` quantifies that geometry and tracks
empirical niche boundary exceedance across space, time and climate models.

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

`climniche` reports four quantities from one niche relative decomposition.
Field names in fitted R objects use snake_case; figures and reports use the
formal names below.

- Climatic Displacement (`climate_change_amount`): distance between current
  and future conditions in the fitted sensitivity weighted climatic space.
- Niche Distance Shift (`niche_distance_change`): signed change in distance
  from the current realised climatic niche centre.
- Climatic Reconfiguration (`climate_reconfiguration`): non radial component
  of climatic displacement. It is derived from Climatic Displacement and Niche
  Distance Shift rather than estimated independently.
- Niche Boundary Exceedance (`niche_boundary_exceedance`): positive excess
  beyond the chosen weighted quantile of current reference cell distances from
  the realised niche centre.

Niche Boundary Exceedance is the boundary outcome. The other three quantities
describe the magnitude and geometry of the climatic displacement that produces
it.

## Basic use

```r
library(climniche)

sim <- simulate_climniche()

fit <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  sensitivity = sim$sensitivity
)

climniche_summary(fit)
climniche_report(fit, species = "example species")
plot_climniche_summary_figure(fit)
```

`fit_climniche_raster()` and `fit_climniche_terra()` accept binary reference
rasters and continuous SDM suitability rasters. Continuous values remain
weights; `occupied_threshold` only sets values at or below the cutoff to zero.
`domain` limits the cells evaluated, while `study_region` adds an optional
boundary to maps.

## Through time

`fit_climniche_series()` holds the fitted current niche centre, climatic
weighting matrix and empirical boundary fixed across future periods and climate
models. For spatial series, future missing cells do not alter this reference;
comparisons use cells available in every projection.

```r
future <- lapply(c(0.25, 0.50, 0.75, 1), function(fraction) {
  sim$current + fraction * (sim$future_away - sim$current)
})

series <- fit_climniche_series(
  current = sim$current,
  future = future,
  time = c(2030, 2050, 2070, 2090),
  occupied = sim$occupied,
  sensitivity = sim$sensitivity
)

climniche_range_summary(series)
climniche_departure(series)
plot_climniche_time(series)
```

Range summaries separate the weighted fraction beyond the empirical niche
boundary from the mean relative exceedance among those cells. Their product is
the weighted mean relative exceedance across the analysed range. This is a
range level summary of Niche Boundary Exceedance rather than an additional
cell level metric.
Optional aggregation and raster cell area weights remain separate from the
reference weights used to estimate the realised niche.

## Worked example

The [Examples](https://bohao0813.github.io/climniche/articles/climniche-examples.html)
page applies the four metrics to European anchovy in the Mediterranean Sea
using OBIS occurrences, Bio-ORACLE v3 climate layers and continuous SDM
suitability weights under SSP2-4.5.

## Contributor

Contributions are welcome. Please make a pull request or contact
bohao.he@polimi.it.
