context("Test utility functions")

test_that("Making unique names", {
    dupNames = c('a-b', 'c-d', 'a-b')
    newNames = makeUniqueNames(dupNames)
    expect_equal(newNames, c('a-b.1', 'c-d', 'a-b.2'))
})
