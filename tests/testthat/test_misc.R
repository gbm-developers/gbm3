context("miscellaneous tests")

test_that("predicts correctly on unknown levels (issue #18)", {
    d <- data.frame(x=as.factor(1:20), y=1:20)

    train <- d[1:10,]
    test  <- d[11:20,]

    ## try this 10 times
    ## this isn't a deterministic bug
    for (ind in seq_len(10)) {
        g <- gbm(y ~ x,
                 distribution="gaussian",
                 bag.fraction=1,
                 data=train, 
                 n.trees=1,
                 shrinkage=1,
                 n.minobsinnode=1)
        p <- predict(g, newdata=test, n.trees=1) - g$initF
        if (sum(abs(p)) != 0) break
    }

    expect_equal(sum(abs(p)), 0)
})
