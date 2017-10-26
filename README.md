gbm3: generalized boosted models
----0------------------------

[![Build Status](https://travis-ci.org/gbm-developers/gbm.svg?branch=master)](https://travis-ci.org/gbm-developers/gbm3)
[![Coverage Status](https://coveralls.io/repos/gbm-developers/gbm/badge.svg?branch=master&service=github)](https://coveralls.io/github/gbm-developers/gbm3?branch=master)

Originally written by Greg Ridgeway, added to by various authors,
currently maintained by Harry Southworth.  Development is discussed
--- somewhat --- at https://groups.google.com/forum/#!forum/gbm-dev .

This is the shiny new gbm3 package that is not backwards compatible, but 
is fast and parallel and --- to some extent --- developed.

Non-production releases (bug fixes, mostly) will be released via the GitHub
release workflow. To install from GitHub, first install `devtools` from CRAN:

```R
install.packages("devtools")
```

Then install `gbm3` from GitHub:

```R
library("devtools")
install_github("gbm-developers/gbm3")
```
