# Mediterranean anchovy case study

`mediterranean_anchovy_case_study.R` generates the fitted case study and its
main figures. `mediterranean_anchovy_priority.R` compares larger positive Niche
Distance Shift with lower Climatic Displacement.
`mediterranean_anchovy_contributions.R` maps the fitted climatic
contributions. `mediterranean_anchovy_time_series.R` evaluates SSP2-4.5 at
2030, 2050, 2070 and 2090 against one fixed current niche reference. Run the
scripts from the package root.

The main case-study script downloads OBIS occurrence records and Bio-ORACLE v3
environmental layers. The time-series script reuses its fitted suitability
surface. `mediterranean_iho_mrgid1905.geojson` is the IHO Mediterranean Sea
area (Marine Regions MRGID 1905) used by all four scripts.

Bio-ORACLE v3 should be cited as Assis et al. (2024), *Global Ecology and
Biogeography* 33, e13813
([doi:10.1111/geb.13813](https://doi.org/10.1111/geb.13813)). OBIS records
should be cited with their source dataset metadata and the
[OBIS database](https://obis.org/).

The boundary is derived from Marine Regions, *IHO Sea Areas, version 3*
([doi:10.14284/323](https://doi.org/10.14284/323)). It was repaired, projected
to EPSG:3035, simplified with a 500 m tolerance and returned to WGS 84. The
simplified geometry and the source geometry selected the same cells at the
example raster resolution. SHA-256:
`21FFA785F1D3CDC20CBCBF398B46C42C0F4923306A3DF20E6910CB4468D63479`.

Raw downloads stay under `../../data-raw/`. Model outputs, including any
manuscript PDF/SVG/TIFF copies, are written under `output/`. Neither location is
included in the CRAN source package. The prepared tables used by the vignette,
including the fitted climatic metric weights, are stored under
`inst/extdata/mediterranean_anchovy/`.

Set `CLIMNICHE_CASE_RUN` to reuse a named output directory when running the
downstream contribution, screening and time-series scripts.
