context('UUID generator')

UUID <- sapply(1:1000, function(x) {sampleID()})

test_that('IDs have the right format', {
    expect_true(all(grepl('[a-f0-9]{8}-[a-f0-9]{4}-4[a-f0-9]{3}-[89aAbB][a-f0-9]{3}-[a-f0-9]{12}', UUID)))
})

test_that('IDs are unique', {
    expect_equal(length(UUID), length(unique(UUID)))
})