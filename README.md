gbm: gradient boosted models
----------------------------

[![Build Status](https://travis-ci.org/gbm-developers/gbm.svg?branch=master)](https://travis-ci.org/gbm-developers/gbm)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/gbm)](https://cran.r-project.org/web/packages/gbm/)

Originally written by Greg Ridgeway, added to by various authors,
currently maintained by Harry Southworth.  Development is discussed
--- somewhat --- at https://groups.google.com/forum/#!forum/gbm-dev .

Non-production releases (bug fixes, mostly) will be released via the GitHub
release workflow. To install from GitHub, first install `devtools` from CRAN:

```R
install.packages("devtools")
```

Then install `gbm` from GitHub:

```R
library("devtools")
install_github("gbm-developers/gbm")
```

Note that we are currently in the middle of a fairly major tidy up, so
you should worry a bit about this code here...
