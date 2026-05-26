<p align="center" style="margin-top: 28px; margin-bottom: 24px;">
  <img src="man/figures/climniche-hex.svg" width="190" alt="climniche logo" />
</p>

# climniche

[![R-CMD-check](https://github.com/Bohao0813/climniche/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Bohao0813/climniche/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Bohao0813/climniche/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Bohao0813/climniche/actions/workflows/pkgdown.yaml)

Website: <https://bohao0813.github.io/climniche/>

Pages: [Examples](https://bohao0813.github.io/climniche/articles/climniche-examples.html) |
[Reference](https://bohao0813.github.io/climniche/reference/index.html) |
[News](https://bohao0813.github.io/climniche/news/index.html)

![Niche climate exposure concept](man/figures/niche-climate-exposure.svg)

`climniche` quantifies niche climate exposure: the amount and direction of
projected climate change at cells currently associated with a taxon, measured
relative to its realised climatic niche.

`climniche` separates four quantities in exposure maps:

- `climate_change_amount`: how far climate moves between the current and future
  period;
- `niche_distance_change`: whether future climate becomes closer to or farther
  from the realised niche centre;
- `composition_change`: change in climate composition that does not mainly
  alter distance from the niche centre;
- `outside_niche_exceedance`: how far future climate lies beyond an empirical
  boundary of current niche conditions.

Together, these quantities distinguish local environmental change from change
that moves current occurrence, range, or SDM cells outside the realised niche.

## Inputs and outputs

`climniche` accepts extracted environmental matrices and aligned `raster` or
`terra` layers. The current reference cells can be supplied as occurrences,
a range raster, or a continuous SDM suitability raster with a chosen threshold.

The outputs are cell-level tables, maps, exposure classes, variable
contributions, report text, and figure data.

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
plot_climniche_diagram(fit)
```

For spatial data, use `fit_climniche_raster()` or `fit_climniche_terra()`.
Both functions accept binary and continuous occupied rasters. `domain` can be
used to restrict the calculation to a study area such as a marine region or a
modelled accessible area.
