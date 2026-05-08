## Resubmission

This is a resubmission of gbm3 3.0.1. The package was previously archived on
2024-02-06 after check problems were not corrected in time.

In this resubmission I have:

* updated the package version and Date field;
* updated stale URLs in the documentation;
* removed stale Travis CI metadata from the CRAN build;
* added GitHub Actions checks for Linux, macOS, and Windows;
* cleaned local build artifacts from the source package;
* reduced slow tests while preserving the tested behavior.

## Test environments

* Local Windows 11 x64, R 4.5.2:
  `R CMD check --as-cran --no-manual`
* macOS arm64, R 4.6.0 patched, macbuilder:
  OK
* GitHub Actions:
  TODO: add results for macOS, Windows, Ubuntu R-devel, Ubuntu release,
  and Ubuntu oldrel after the workflow has run.
* Optional external checks:
  TODO: add R-hub/win-builder results if run before submission.

## R CMD check results

Local Windows 11 x64, R 4.5.2:

0 errors | 0 warnings | 3 notes

The notes were:

* New submission / package was archived on CRAN.
* README.md or NEWS.md could not be checked because pandoc was not detected by
  the local check process.
* One example exceeded 5 seconds of CPU/elapsed time.

## Previous CRAN archive issue

The archived CRAN check results for gbm3 3.0 on 2024-02-06 reported an ERROR on
`r-devel-linux-x86_64-debian-clang` because required/suggested packages were not
available for checking in that environment. The package now declares its required
and suggested dependencies in DESCRIPTION, including Rcpp, testthat, rmarkdown,
knitr, and MASS.
