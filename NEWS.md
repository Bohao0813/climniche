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
