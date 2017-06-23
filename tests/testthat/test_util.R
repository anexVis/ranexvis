context("Test utility functions")

test_that("Making unique names", {
    dupNames = c('a-b', 'c-d', 'a-b')
    newNames = makeUniqueNames(dupNames)
    expect_equal(newNames, c('a-b.1', 'c-d', 'a-b.2'))
})

test_that("Removing Ensembl version", {
    withVersion = c("ENSG00000002587.5", "ENSG00000125430.4", "ENSG00000176022.3", "ENSG00000027847.9")
    withoutVersion = c("ENSG00000002587", "ENSG00000125430", "ENSG00000176022", "ENSG00000027847")
    expect_equal(removeEnsemblVersion(withVersion), withoutVersion)
})

test_that("HGNC --> Ensembl IDs", {
    geneSymbol = c('HS3ST1', 'HS3ST3B1', 'B3GALT6', 'B4GALT7')
    geneEnsembl = c("ENSG00000002587.5", "ENSG00000125430.4", "ENSG00000176022.3", "ENSG00000027847.9")
    expect_equal(hgnc2ensembl(geneSymbol), removeEnsemblVersion(geneEnsembl))
})

test_that("Ensembl gene IDs --> HGNC (local query)", {
    geneSymbol = c('HS3ST1', 'HS3ST3B1', 'B3GALT6', 'B4GALT7')
    geneEnsembl = c("ENSG00000002587", "ENSG00000125430", "ENSG00000176022", "ENSG00000027847")
    expect_equal(ensembl2hgnc(geneEnsembl), geneSymbol)
})

## Disable remote query
# test_that("Ensembl gene IDs --> HGNC (remote query)", {
#     geneSymbol = c('HS3ST1', 'HS3ST3B1', 'B3GALT6', 'B4GALT7')
#     geneEnsembl = c("ENSG00000002587", "ENSG00000125430", "ENSG00000176022", "ENSG00000027847")
#     expect_equal(ensembl2hgnc.remote(geneEnsembl), geneSymbol)
# })
