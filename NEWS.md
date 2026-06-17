# climniche 0.2.0

* Reframed user-facing documentation around a niche-relative decomposition of
  projected climate exposure rather than four independent metrics.
* Reworked the Mediterranean anchovy vignette to start from prepared climniche
  inputs and to separate Bio-ORACLE layer metadata, presence-background SDM
  metadata, predictor screening metadata and climniche fit settings.
* Added AUC to the example SDM metadata and used equal numbers of presence and
  background cells in the prepared example.
* Updated figures, reports and package-site examples to use the formal names
  Climatic Displacement, Niche Distance Shift, Climatic Reconfiguration and
  Niche Boundary Exceedance for the reported quantities.
* Removed the Mediterranean anchovy diagram from the vignette because the
  two-dimensional projection did not provide a clear enough explanation of the
  niche-relative decomposition.
* Added publication-oriented PDF outputs for the Mediterranean anchovy map and
  summary figures used in the vignette.
* Refreshed the Mediterranean anchovy example with prepared Bio-ORACLE v3
  marine climate rasters and continuous suitability weights.
* Reworked the summary figure to use weighted exposure-plane bins, weighted
  quantity distributions, normalised variable contribution shares and all
  five possible hierarchical exposure classes, including zero proportions.
* Restored the package default classification settings in the anchovy example;
  example thresholds are no longer adjusted to produce a chosen number of
  classes.
* Added `plot_climniche_maps()` and extended map controls for colour limits,
  palettes, legends, plotting extent and longitude-latitude labels.
* Revised the Mediterranean maps to use a dissolved sea boundary, a distinct
  land mask, aligned geographic axes and horizontal legends.
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
* Classification thresholds can now be set by the user through `tolerance`,
  `stable_climate_change`, `stable_reconfiguration`,
  `boundary_exceedance_tolerance` and the associated quantile arguments.
* Fitted objects now store the effective values used for classification in
  `classification_settings`.

# climniche 0.0.1

* Initial CRAN release.
* Provided matrix, raster and terra workflows for niche climate exposure.
* Added report, summary and figure functions for cell-level exposure metrics.
