## Resubmission

This is a resubmission of gbm3 3.0.1. The package was previously archived on
2024-02-06 after check problems were not corrected in time.

In this resubmission I have:

* updated the package version and Date field;
* updated stale URLs in the documentation;
* removed stale Travis CI metadata from the CRAN build;
* added GitHub Actions checks for Linux, macOS, and Windows;
* cleaned local build artifacts from the source package;
* removed undefined behavior detected by sanitizer checks;
* reduced slow tests while preserving the tested behavior;
* skipped long-running integration-style tests during valgrind checks to avoid
  CRAN's special-check timeout. These tests continue to run during ordinary
  checks.

## Test environments

* Local Windows 11 x64, R 4.6.0:
  `R CMD check --as-cran --no-manual`
* macOS arm64, R 4.6.0 patched, macbuilder:
  OK
* R-hub:
  valgrind, clang-asan, clang-ubsan, gcc-asan, ubuntu-clang, ubuntu-gcc12,
  windows, clang20, gcc15: OK
* CRAN incoming pretest special checks:
  clang-san and gcc-san: OK. valgrind initially timed out in long-running
  integration-style tests with no Valgrind memory errors reported; those tests
  are now skipped only during Valgrind checks.
* GitHub Actions R-CMD-check:
  macos-latest, ubuntu-latest (devel), ubuntu-latest (oldrel-1),
  ubuntu-latest (release), windows-latest (release): OK


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
