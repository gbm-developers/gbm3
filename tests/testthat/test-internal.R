context("testing miscellaneous helper functions")

test_that("warnNoVariation warns if a variable has no variation", {
    x <- c(1.2, 1.2)
    
    expect_warning(warnNoVariation(x, 1, 'test'),
                   "variable 1: test has no variation",
                   fixed=TRUE)
})

test_that("warnNoVariation passes OK if variable does vary", {
    x <- c(0, 1)

    expect_warning(warnNoVariation(x, 1, 'test'), regexp=NA)
})

test_that("warnNoVariation passes OK if variable does vary (with NA)", {
    x <- c(0, 1, NA)

    expect_warning(warnNoVariation(x, 1, 'test'), regexp=NA)
})
