## Test environments

* local Windows, R release

## R CMD check results

0 errors | 0 warnings | 1 note

* Non-staged installation was used for the local Windows check because the
  sandboxed filesystem prevented R from renaming the staged installation
  directory. The package installed, loaded, ran its tests and rebuilt its
  vignettes successfully under this check.

This update follows the initial CRAN release. It revises terminology used in
user-facing output, adds user-set classification thresholds and treats
continuous SDM suitability as reference weights.

## Downstream dependencies

There are no known downstream dependencies.
