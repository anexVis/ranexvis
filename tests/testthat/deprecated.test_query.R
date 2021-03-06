library(futile.logger)
context("Setting up may take some time. Data retrieval")
flog.appender(appender.file("test_query.log"), name="log")

# Need to set up before querying

# setup()

test_that("Container was loaded with data", {
    expect_true(exists("geneList", container))
    expect_true(exists("sampleMetadata", container))
    expect_true(exists("/toil-rsem/gene/tpm/expressionMatrix", container))
})

test_that("Database path is set correctly.", {
    expect_equal(dbpath[['gtex']], '/var/www/html/data/GTEx_V6.h5')
})

test_that("Retrieval of gene list.", {
    # expect_warning(genelist = getGeneList(expect='df'),'HDF5 on bit64conversion')
    # expect_warning(jsonlist = getGeneList(db='gtex', expect='json'), 'HDF5 on bit64conversion')
    genelist = getGeneList("gtex", c("EnsemblID", "HGNC"), expect="df", read.from.redis=FALSE)
    expect_equal(names(genelist) , c("EnsemblID", "HGNC"))
})

test_that("Retrieval of gene list without Ensembl version number", {
    genelist = getGeneList("gtex", c("EnsemblID", "HGNC"), withEnsemblVersion=FALSE,expect="df", read.from.redis = FALSE)
    expect_equal(names(genelist) , c("EnsemblID", "HGNC"))
})

test_that("Retrieval of gene sets", {
    genesets = getGeneSets(db="gtex", processing="toil-rsem", expect="json",read.from.redis=FALSE)
    flog.info("Gene sets JSON: ", jsonlite::prettify(genesets), name='log', capture=TRUE)
})

test_that("Retrieval of sample grouping list.", {
   sampleGroupingList = getSampleGroupingList(expect='df',read.from.redis=FALSE)
})

test_that("Retrieval of sample metadata.", {
   sampleMetadata = getSampleMetadata("gtex", cols=c("SAMPID", "SMTS"), expect="dt",read.from.redis=FALSE)
   expect_equal(ncol(sampleMetadata), 2)
})

test_that("Retrieval of sample metadata by group", {
    sampleMetadata = getSampleMetadataByGroup("Bladder", sampleGrouping="SMTS", cols=c("SAMPID", "SMTS","SMUBRID"), expect="dt", read.from.redis=FALSE)
    flog.info("Sample metadata:", sampleMetadata[1:10,], name='log', capture=TRUE)
    expect_equal(ncol(sampleMetadata),3)
    expect_equal(sampleMetadata$SMTS,rep("Bladder", nrow(sampleMetadata)))
})
# obsolete test
# test_that("Reading 1D array from HDF5", {
#    db = 'gtex'
#    geneList = readCharacterArray(dbpath[[db]], "/genes/EnsemblID")
#    expect_equal(class(geneList), "character")
#    geneListDf= readCharacterArray(dbpath[[db]], "/genes/EnsemblID", colname="EnsemblID")
#    expect_equal(class(geneListDf), "data.frame")
#    expect_equal(nrow(geneListDf), length(geneList))
#    # expect_equal(nrow(geneList) , 60498)
# })
#
test_that("Retrieval of gene expression matrix.", {
   # selectedGenes = c("HS3ST1", "HS3ST3B1")
   # selectedEnsembls = c("ENSG00000002587.5", "ENSG00000125430.4") # HS3ST1, HS3STB1
   selectedEnsembls = c("ENSG00000176022.3", "ENSG00000027847.9") # B3GALT6, B4GALT7
   # expect_error()   # when db/processing/unit are not found

   expr1 = getGeneExpressionMatrix(genes = selectedEnsembls,
                                  sampleGroups="Brain",
                                  sampleGrouping="SMTS",
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable',
                                  read.from.redis=FALSE
                                  )
   expect_false(is.null(expr1))
   expect_equal(ncol(expr1),2)

   ### Note that the returned matrix will have EnsemblID stripped off version number
   expect_equal(colnames(expr1),removeEnsemblVersion(selectedEnsembls))
   flog.info("Expression of 2 genes in the first 10 Brain samples:", expr1[1:10,], name='log', capture=TRUE)


   expr2 = getGeneExpressionMatrix(genes = selectedEnsembls,
                                  sampleGroups=c("Brain", "Liver"),
                                  sampleGrouping="SMTS",
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable',
                                  read.from.redis=FALSE
                                  )
   expect_equal(ncol(expr2),2)
   expect_equal(colnames(expr1), colnames(expr2))
   expect_true(nrow(expr1) < nrow(expr2))
   expect_equal(nrow(expr1), 1146)  # 1146 brain samples
   expect_equal(cor(expr1,method="pearson")[1,1], 1)
   expect_true(cor(expr1,method="pearson")[1,2] > 0.5) # these two genes are known to be highly correlated

   # querying non-existing genes
   expect_message(getGeneExpressionMatrix(genes = "aliengene",
                                  sampleGroups=c("Brain", "Liver"),
                                  sampleGrouping="SMTS",
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable',
                                  read.from.redis=FALSE
                                  ), "No matching record found.")
})

test_that("Retrieval of gene expression matrix in JSON string.", {
    selectedEnsembls = c("ENSG00000176022.3", "ENSG00000027847.9") # B3GALT6, B4GALT7
    expr_json = getGeneExpressionMatrix(genes = selectedEnsembls,
                                  sampleGroups="Bladder",
                                  sampleGrouping="SMTS",
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='json',
                                  read.from.redis=FALSE
                                  )
    flog.info("Expression of 2 genes in Bladder samples, JSON format:", jsonlite::prettify(expr_json), name='log', capture=TRUE)
})

test_that("Retrieval of scatter data in JSON string.", {
    # B3GALT6, B4GALT7
    expr_json = getScatterData(x = "ENSG00000176022.3",
                               y = "ENSG00000027847.9",
                               sampleGroups="Bladder",
                               sampleGrouping="SMTS",
                               db = "gtex",
                               processing = "toil-rsem",
                               unit="tpm",
                               read.from.redis=FALSE
                               )
    flog.info("Scatter data in JSON format:", jsonlite::prettify(expr_json), name='log', capture=TRUE)
})

test_that("Retrieval of scatter data for identical pair.", {
    # B3GALT6, B4GALT7
    expr_json = getScatterData(x = "ENSG00000176022.3",
                               y = "ENSG00000176022.3",
                               sampleGroups="Bladder",
                               sampleGrouping="SMTSD",
                               db = "gtex",
                               processing = "toil-rsem",
                               unit="tpm",
                               read.from.redis=FALSE
                               )
    flog.info("Scatter data for identical pair:", jsonlite::prettify(expr_json), name='log', capture=TRUE)
})
