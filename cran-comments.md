## Test environments

* local Windows 11, R 4.6.0 (`R CMD check --as-cran --no-manual`)
* GitHub Actions, Ubuntu, R release

## R CMD check results

0 errors | 0 warnings | 0 notes

Version 0.3.8 is the first CRAN update since 0.0.1. It adds continuous
reference weights, configurable niche descriptors, climate variable
preprocessing, raster and terra workflows, fixed reference projections,
ordered time series, climatic-variable contributions, ecological screening,
and report and figure methods. The original field names remain available as
aliases.

The package license has changed from GPL (>= 3) to MIT. The copyright holder
is unchanged. Third-party sources used to prepare the case-study materials are
identified in `inst/COPYRIGHTS`; unmodified OBIS occurrence records and Marine
Regions geometry are not distributed with the package.

## Downstream dependencies

There are no known downstream dependencies.
