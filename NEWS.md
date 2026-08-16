# climniche 0.3.8

* Separated first observed Niche Boundary Exceedance from persistent
  exceedance onset. A first exceedance can now occur in the final supplied
  projection without being treated as missing.
* Added primary fields for the fraction of supplied projections beyond the
  niche boundary while retaining the existing temporal field names.
* Revised the temporal example to map first boundary exceedance across all
  supplied projection years.

# climniche 0.3.7

* Clarified the mathematical relation among Climatic Displacement, Niche
  Distance Shift and the derived Climatic Reconfiguration term, and separated
  that geometry from the future-state Niche Boundary Exceedance statistic.
* Corrected the upper endpoint of weighted niche percentiles and projected
  small numerical negative eigenvalues to the nearest positive semidefinite
  climatic metric.
* Fitted spatial current references before applying future missing-value masks
  and preserved named reference weights when projecting subsets of cells.
* Required fitted CENFA components to use their original variables and
  standardisation.
* Made `summary.climniche_fit()` reference-weighted by default and aligned the
  default Pareto filter with the selected exposure direction.
* Added executable controlled-geometry and Mediterranean input examples, plus
  a documented IHO Mediterranean boundary for rebuilding the case study.

# climniche 0.3.6

* Added a Bio-ORACLE SSP2-4.5 European anchovy time-series example with
  range summaries, persistent Niche Boundary Exceedance maps and the
  series-report workflow.
* Removed an unused fill scale from single-model time plots.
* Standardised user-facing names for range, temporal and climatic contribution
  summaries without removing existing R object fields.
* Added `squared_niche_distance_change` as the primary climatic contribution
  field while retaining `niche_potential_change` as an exact compatibility
  alias.

# climniche 0.3.5

* Reframed Pareto output as ecological screening of climate exposure.
  `exposure_direction` distinguishes larger positive Niche Distance Shift from
  lower Climatic Displacement.
* Added `pareto_depth_score` as the primary Pareto depth field while retaining
  `relative_priority` as an exact compatibility alias.
* Reworked the Mediterranean anchovy example to compare positive Niche Distance
  Shift with low Climatic Displacement in paired decision planes and maps.

# climniche 0.3.4

* Added two-objective Pareto depth for climate exposure screening, including
  decision planes and spatial Pareto maps.
* Added a European anchovy screening example and clarified that summary figure
  panel (c) decomposes squared niche distance change rather than SDM variable
  importance.
* Added cell-level dominant climatic contributions, dominance shares and
  spatial maps, with a Mediterranean anchovy example.
* Added complete variable legends for contribution maps and aligned the
  Mediterranean example to one six-variable climate set.
* Added Pareto summary diagnostics, including first-front size and objective
  rank correlation.

# climniche 0.3.3

* Added `fit_climniche_reference()` and `project_climniche()` to estimate the
  current realised climatic niche once and reuse it across projections.
* Added `fit_climniche_series()` for ordered future periods, climate models and
  scenarios. Series share one effective descriptor tolerance. Spatial
  references are fitted independently of future missing cells, while
  projections use common finite cell support.
* Added weighted range summaries of Niche Boundary Exceedance: its weighted
  fraction, conditional relative magnitude and range mean. Reference,
  aggregation and raster cell area weights remain distinct.
* Added continuous summaries of persistent exceedance onset, time weighted and
  cumulative exceedance, return below the boundary and maximum interval
  increase.
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
* Added executable examples, a table of case-study climatic metric weights and
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

* Focused fitted objects, reports and figures on the four continuous reported
  quantities and the two descriptor fields.
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
  fitted niche centre, empirical radial boundary, niche percentiles and
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
