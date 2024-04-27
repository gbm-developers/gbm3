gbm3: generalized boosted models
----0------------------------

<!-- [![Build Status](https://travis-ci.org/gbm-developers/gbm.svg?branch=master)](https://travis-ci.org/gbm-developers/gbm3) -->
<!--  [![Coverage Status](https://coveralls.io/repos/gbm-developers/gbm/badge.svg?branch=master&service=github)](https://coveralls.io/github/gbm-developers/gbm3?branch=master) -->

Originally written by Greg Ridgeway between 1999-2003, added to by various 
authors, extensively updated and polished by James Hickey in 2016, survival 
models greatly improved by Terry Therneau in 2016, and currently 
maintained by Greg Ridgeway.
Development is discussed --- somewhat --- at
https://groups.google.com/forum/#!forum/gbm-dev .

This is the shiny new gbm3 package that is not backwards compatible with R code
calling the original gbm package, but is fast and parallel and developed.

Non-production releases (bug fixes, mostly) will be released via the GitHub
release workflow. To install from GitHub, first install `remotes` from CRAN:

```R
install.packages("remotes")
```

Then install `gbm3` from GitHub:

```R
remotes::install_github("gbm-developers/gbm3")

# or to ensure your got everything
remotes::install_github("gbm-developers/gbm3", build_vignettes = TRUE, force = TRUE)
```
