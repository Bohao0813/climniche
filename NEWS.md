# climniche 0.4.0

* Focused fitted objects, reports and figures on the four continuous exposure
  metrics and the two descriptor fields.
* Updated the Mediterranean anchovy maps so metric panels are drawn only within
  current suitable habitat while retaining the Mediterranean Sea as the
  analysis domain.
* Revised the summary figure labels, panel tags and metric distributions for
  manuscript use.
* Added PDF outputs for the prepared Mediterranean maps and summary figure.

# climniche 0.3.0

* Revised the Mediterranean anchovy example so the analysis domain is the
  Mediterranean Sea mask, while the continuous SDM suitability raster is used
  only as the current reference weight layer.
* Added `plot_climniche_summary_figure()` and
  `climniche_summary_figure_data()`.
* Updated map and summary figure outputs in PNG and PDF format for the prepared
  Mediterranean example.
* Cleaned README, vignette, report and pkgdown text to use formal metric names
  and more direct methodological wording.
* Added tests for formal plot/report labels, domain handling with continuous
  reference weights and threshold values equal to `occupied_threshold`.

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
