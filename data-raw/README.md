# Mediterranean anchovy case study

`mediterranean_anchovy_case_study.R` generates the fitted case study and its
main figures. `mediterranean_anchovy_priority.R` compares ecological exposure
concern with climatic persistence opportunity.
`mediterranean_anchovy_contributions.R` maps the fitted climatic
contributions. `mediterranean_anchovy_time_series.R` evaluates SSP2-4.5 at
2030, 2050, 2070 and 2090 against one fixed current niche reference. Run the
scripts from the package root.

The main case-study script downloads OBIS occurrence records and Bio-ORACLE v3
environmental layers. The time-series script reuses its fitted suitability
surface. Both expect the IHO Mediterranean Sea boundary (Marine Regions MRGID
1905) outside the package tree at:

```text
../../data-raw/marine_regions/mediterranean_iho_mrgid1905.gpkg
```

Raw downloads stay under `../../data-raw/`. Model outputs, including any
manuscript PDF/SVG/TIFF copies, are written under `output/`. Neither location is
included in the CRAN source package. The prepared tables used by the vignette,
including the fitted sensitivity weights, are stored under
`inst/extdata/mediterranean_anchovy/`.
