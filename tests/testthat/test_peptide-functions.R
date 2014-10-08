context('Peptide functions')

test_that('mass calculation works', {
    expect_error(pepMass('ABCD'))
    expect_equal(pepMass('STY'), c(STY=369.4))
    expect_equal(pepMass('STY', mono=TRUE), c(STY=369.153599))
    expect_equal(pepMass('STY', mono=TRUE, neutral=TRUE), c(STY=351.143039))
    expect_equal(pepMass('STY', modifications=c(10, -5)), c(STY=374.4))
})

test_that('fragmentation pattern works', {
    pattern <- fragPattern('STY')
    expect_equal(names(pattern), c('ion', 'index', 'mz'))
    expect_equal(sum(pattern$mz), 3728.675408)
    expect_equal(unique(pattern$index), 1:2)
    expect_true(all(pattern$ion %in% c("a", "a\u00BA", "b", "b\u00BA", "c", "c\u00BA", "x", "x\u00BA", "y", "y\u00BA", "z", "z\u00BA")))
    
    pattern <- fragPattern('STY', list(NULL, 10, NULL), ions='aby', neutralLosses=FALSE)
    expect_equal(sum(pattern$mz), 993.476026)
    expect_true(all(pattern$ion %in% c("a", "b", "y")))
})