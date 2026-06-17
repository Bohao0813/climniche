# Mediterranean anchovy case study

`mediterranean_anchovy_case_study.R` generates the prepared CSV files and the
PNG, PDF, SVG and TIFF figures used by the package vignette. Run it from the
package root.

The script downloads OBIS occurrence records and Bio-ORACLE v3 environmental
layers. It expects the IHO Mediterranean Sea boundary (Marine Regions MRGID
1905) at:

```text
data-raw/marine_regions/mediterranean_iho_mrgid1905.gpkg
```

Raw downloads and model outputs remain under `data-raw/` and `output/`; neither
directory is included in the CRAN source package. The prepared, compact tables
used by the vignette are stored under `inst/extdata/mediterranean_anchovy/`.
