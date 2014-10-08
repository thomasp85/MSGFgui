context('Data tables')

test_that('table gets loaded', {
    expect_equal(dim(getAAtable()), c(20, 15))
    expect_equal(dim(getAdductTable()), c(32, 4))
})