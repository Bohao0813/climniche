# climniche 0.3.4

* Added two-objective Pareto depth for climate exposure screening, including
  decision planes and spatial priority maps.
* Added a European anchovy priority example and clarified that summary figure
  panel (c) decomposes change in niche potential rather than SDM variable
  importance.
* Added cell-level dominant climatic contributions, dominance shares and
  spatial maps, with a Mediterranean anchovy example.
* Added Pareto summary diagnostics, including first-front size and objective
  rank correlation.

# climniche 0.3.3

* Added `fit_climniche_reference()` and `project_climniche()` to estimate the
  current realised climatic niche once and reuse it across projections.
* Added `fit_climniche_series()` for ordered future periods, climate models and
  scenarios. Series share one effective descriptor tolerance. Spatial
  references are fitted independently of future missing cells, while
  projections use common finite cell support.
* Added a range level decomposition of Niche Boundary Exceedance into exposed
  fraction, conditional relative exceedance and their range mean. Reference,
  aggregation and raster cell area weights remain distinct.
* Added continuous summaries of departure onset, duration, cumulative
  exceedance, re-entry and maximum interval increase.
* Added climate model agreement for ensembles with at least two models, time
  series reports and dynamic time and map figures.
* Documented the angular identity underlying Climatic Reconfiguration.

# climniche 0.3.2

* Matched named climate rows, variables, reference weights, centres,
  sensitivities, metric matrices and compatible CNFA components before fitting.
* Made correlation filtering deterministic and used one inverse weighted
  empirical distribution definition for binary and continuous reference
  weights.
* Added stricter validation for fitted climate spaces and factor metrics. A
  requested factor metric now stops when its CNFA components are incomplete.
* Corrected weighted summaries and plots, single-layer raster handling,
  non-finite raster masks, spatial layer names and original raster cell indices.
* Added executable examples, a table of case-study sensitivity weights and
  extended tests for alignment, weighting and invalid inputs.
* Kept manuscript PDF figures in the repository while excluding them from the
  source package built for CRAN.

# climniche 0.3.1

* Added default climate-variable preprocessing in `fit_climniche()` and the
  raster workflows. The preprocessing removes near-zero variance variables and
  highly correlated variables before metric fitting.
* Added preprocessing metadata to fitted objects.
* Documented compatibility with CENFA-style `cnfa` objects without requiring
  the CENFA package at installation time.
* Refreshed the Mediterranean anchovy example outputs under the current
  analysis workflow.
* Revised the `fit_climniche()` formula documentation for clearer pkgdown
  rendering.

# climniche 0.3.0

* Focused fitted objects, reports and figures on the four continuous exposure
  metrics and the two descriptor fields.
* Added `plot_climniche_summary_figure()` and
  `climniche_summary_figure_data()`.
* Added map controls for binary and continuous SDM reference rasters, study
  region boundaries and longitude-latitude labels.
* Updated report, summary, map and vignette outputs to use formal metric names.
* Added tests for descriptor thresholds, continuous reference weights and
  raster domain handling.

# climniche 0.2.0

* Added the Mediterranean European anchovy example with Bio-ORACLE v3 marine
  climate layers, prepared future projections and continuous SDM suitability
  weights.
* Added `radial_direction` and `boundary_status` descriptors to fitted objects,
  cell tables and raster workflows.
* Added `plot_climniche_maps()` and map controls for colour limits, palettes,
  legends, plotting extent, study-region boundaries and longitude-latitude
  labels.
* Changed the package license to MIT.

# climniche 0.1.1

* Continuous SDM suitability inputs are now used as reference weights in the
  fitted niche centre, empirical niche boundary, niche percentiles and
  current-scope summaries.
* `occupied_threshold` removes low reference weights while preserving the
  original continuous values above the threshold.
* Raster and terra workflows now pass continuous suitability rasters through to
  the matrix workflow as weights instead of converting them to binary cell
  indices.

# climniche 0.1.0

* Added the primary R object fields `climate_reconfiguration` and
  `niche_boundary_exceedance`.
* Retained `composition_change` and `outside_niche_exceedance` as compatibility
  aliases for existing code.

# climniche 0.0.1

* Initial CRAN release.
* Provided matrix, raster and terra workflows for niche climate exposure.
* Added report, summary and figure functions for cell-level exposure metrics.
