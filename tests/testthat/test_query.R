test_that("Database path is set correctly.", {
    expect_equal(dbpath[['gtex']], '/var/www/html/data/GTEx_V6.h5')
})

test_that("Retrieval of gene list.", {
    # expect_warning(genelist = getGeneList(expect='df'),'HDF5 on bit64conversion')
    # expect_warning(jsonlist = getGeneList(db='gtex', expect='json'), 'HDF5 on bit64conversion')
    genelist = getGeneList("gtex", c("EnsemblID", "HGNC"), expect="df")
})

test_that("Retrieval of sample grouping list.", {
    sampleGroupingList = getSampleGroupingList(expect='df')
})

test_that("Retrieval of sample metadata.", {
    sampleMetadata = getSampleMetadata("gtex", cols=c("SAMPID", "SMTS"), expect="dt")
    expect_equal(ncol(sampleMetadata), 2)
})

test_that("Reading 1D array from HDF5", {
    db = 'gtex'
    geneList = readCharacterArray(dbpath[[db]], "/genes/EnsemblID")
    expect_equal(class(geneList), "character")
    geneListDf= readCharacterArray(dbpath[[db]], "/genes/EnsemblID", colname="EnsemblID")
    expect_equal(class(geneListDf), "data.frame")
    expect_equal(nrow(geneListDf), length(geneList))
    # expect_equal(nrow(geneList) , 60498)
})
test_that("Retrieval of gene expression matrix.", {
    # selectedGenes = c("HS3ST1", "HS3ST3B1")
    selectedEnsembls = c("ENSG00000002587.5", "ENSG00000125430.4")
    # expect_error()   # when db/processing/unit are not found

    expr = getGeneExpressionMatrix(genes = selectedEnsembls,
                                   sampleGroups=c("Brain"),
                                   sampleGrouping="SMTS",
                                   db = "gtex",
                                   processing = "toil-rsem",
                                   unit="tpm"
                                   )
    expect_equal(nrow(expr),2);
})

