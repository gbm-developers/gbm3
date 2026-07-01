# gbm3: Generalized Boosted Models

<!-- badges: start -->
[![R-CMD-check](https://github.com/gbm-developers/gbm3/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gbm-developers/gbm3/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/gbm-developers/gbm3/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/gbm-developers/gbm3/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

Originally written by Greg Ridgeway between 1999-2003, added to by various 
authors, extensively updated and polished by James Hickey in 2016, survival 
models greatly improved by Terry Therneau in 2016, and currently 
maintained by Greg Ridgeway.
Development is discussed at the
[gbm-dev Google Group](https://groups.google.com/forum/#!forum/gbm-dev).

`gbm3` provides generalized boosted regression models with a newer API than the
original `gbm` package. The package supports regression, classification,
survival models, and learning-to-rank methods, with optional OpenMP
parallelization in the core fitting code.

Documentation and vignettes are available at
<https://gbm-developers.github.io/gbm3/>.

To install the development version from GitHub, first install `remotes`:

```R
install.packages("remotes")
```

Then install `gbm3`:

```R
remotes::install_github("gbm-developers/gbm3")

# or to build vignettes during installation
remotes::install_github("gbm-developers/gbm3", build_vignettes = TRUE, force = TRUE)
```


