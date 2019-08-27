context('Utility functions')

test_that('Standardization works', {
    d <- data.frame(a=c(1,2,3), b=c(2,4,6), c=c(1,1,1))
    ds <- std_normalize(d, c('a','b'))   # c must not be included!
    expect_is(ds, 'data.frame')
    dscmp <- structure(data.frame(a=c(-1,0,1), b=c(-1,0,1), c=c(1,1,1)), means=c(a=2,b=4), stdevs=c(a=1,b=2))
    expect_equal(ds, dscmp)
})
