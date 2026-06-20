<p align="center" style="margin-top: 28px; margin-bottom: 24px;">
  <img src="man/figures/climniche-hex.svg" width="190" alt="climniche logo" />
</p>

# climniche

Climate change becomes exposure when it is read against the climate a species
already occupies. `climniche` measures where each projected climate moves in
that niche: how far, toward or away from the realised niche centre, into a
different climatic combination, or beyond the current boundary.

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
  and future conditions in the fitted sensitivity weighted climatic space.
- Niche Distance Shift (`niche_distance_change`): signed change in distance
  from the current realised climatic niche centre.
- Climatic Reconfiguration (`climate_reconfiguration`): non radial component
  of climatic displacement. It is derived from Climatic Displacement and Niche
  Distance Shift rather than estimated independently.
- Niche Boundary Exceedance (`niche_boundary_exceedance`): positive excess
  beyond the chosen weighted quantile of current reference-cell distances from
  the realised niche centre.

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

For raster workflows, use `fit_climniche_raster()` or `fit_climniche_terra()`.
The reference layer can be binary or continuous. Continuous SDM suitability
values are treated as weights; `occupied_threshold` removes low-suitability
cells but keeps the remaining suitability values continuous. Use `domain` to
restrict the analysis to the accessible area or study region. Map functions can
draw `sf`, `Spatial`, `SpatVector` or x-y boundary objects through
`study_region`.

## Worked example

The [Examples](https://bohao0813.github.io/climniche/articles/climniche-examples.html)
page uses European anchovy in the Mediterranean Sea. OBIS occurrences and
Bio-ORACLE v3 predictors are used to prepare a presence-background SDM. The
SDM cutoff is selected by maximum test-set TSS; cells below the cutoff are
removed from the reference set, and cells above the cutoff keep their
continuous suitability values as reference weights. The example then maps the
four exposure metrics inside current suitable habitat under SSP2-4.5.

## Contributor

Contributions are welcome. Please make a pull request or contact
bohao.he@polimi.it.
