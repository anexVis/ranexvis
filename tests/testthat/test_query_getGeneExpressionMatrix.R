library(futile.logger)
context("Querying gene expression data from redis server")
flog.appender(appender.file("test_query_getGeneExpressionMatrix.log"), name="log")

test_that("Retrieval of gene expression matrix.", {
   # selectedGenes = c("HS3ST1", "HS3ST3B1")
   # selectedEnsembls = c("ENSG00000002587.5", "ENSG00000125430.4") # HS3ST1, HS3STB1
   selectedEnsembls = c("ENSG00000176022.3", "ENSG00000027847.9") # B3GALT6, B4GALT7
   colNames = c("ENSG00000176022", "ENSG00000027847")
   # expect_error()   # when db/processing/unit are not found

   expr1 = getGeneExpressionMatrix(genes = selectedEnsembls,
                                  sampleGroups="Brain",
                                  sampleGrouping="SMTS",
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable'
                                  )
   ### Note that the returned matrix will have EnsemblID stripped off version number
   expect_equal(names(expr1),colNames)

   flog.info("Expression of 2 genes in the first 10 Brain samples:", expr1[1:10,], name='log', capture=TRUE)

   expr2 = getGeneExpressionMatrix(genes = selectedEnsembls,
                                  sampleGroups=c("Brain", "Liver"),
                                  sampleGrouping="SMTS",
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable'
                                  )
   expect_equal(names(expr2),colNames)
   expect_equal(names(expr1), names(expr2))
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
                                  expect='datatable'
                                  ), "No matching record found.")

   expr3 =  getGeneExpressionMatrix(genes = selectedEnsembls,
                                  sampleGroups=c("Brain", "Liver"),
                                  sampleGrouping="SMTS",
                                  sampleMetaFields=c("SMTS", "SMTSD"),
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable'
                                  )
   flog.info("Expression matrix with sample meta data:", expr3[1:10,], name='log', capture=TRUE)
   expect_equal(ncol(expr3),length(selectedEnsembls) + 2)

    expr3j =  getGeneExpressionMatrix(genes = selectedEnsembls,
                                  sampleGroups="Bladder",
                                  sampleGrouping="SMTS",
                                  sampleMetaFields=c("SMTS", "SMTSD"),
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='json'
                                  )
   flog.info("Expression matrix with sample meta data, in JSON:", jsonlite::prettify(expr3j), name='log', capture=TRUE)

})

test_that("Retrieval of expression matrix of single gene", {
   expr1 =  getGeneExpressionMatrix(genes = "ENSG00000176022.3",
                                  sampleGroups=c("Bladder"),
                                  sampleGrouping="SMTS",
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable'
                                  )

   expect_equal(ncol(expr1), 1)
   # flog.info("Expression matrix of single gene:", expr1, name='log', capture=TRUE)
})

test_that("Retrieval of expression matrix of single gene with one-column metadata", {
   expr1 =  getGeneExpressionMatrix(genes = "ENSG00000176022.3",
                                  sampleGroups=c("Bladder"),
                                  sampleGrouping="SMTS",
                                  sampleMetaFields="SMTS",
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable',
                                  read.from.redis = TRUE
                                  )
   expect_equal(names(expr1), c("ENSG00000176022", "SMTS"))

   flog.info("Expression matrix of single gene with metadata:", expr1, name='log', capture=TRUE)
})


test_that("Retrieval of expression matrix of single gene with multiple-column metadata", {
   expr1 =  getGeneExpressionMatrix(genes = "ENSG00000176022.3",
                                  sampleGroups=c("Bladder"),
                                  sampleGrouping="SMTS",
                                  sampleMetaFields=c("SMTS", "SMTSD"),
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable',
                                  read.from.redis = TRUE
                                  )
   expect_equal(ncol(expr1), 3)
   expect_equal(names(expr1), c("ENSG00000176022", "SMTS", "SMTSD"))

   flog.info("Expression matrix of single gene with metadata:", expr1, name='log', capture=TRUE)
})

test_that("Retrieval of non-existent expression matrix", {
    expr0 = getGeneExpressionMatrix(genes = "ENSG0000017",
                                  sampleGroups=c("Bladder"),
                                  sampleGrouping="SMTS",
                                  sampleMetaFields=c("SMTS", "SMTSD"),
                                  db = "gtex",
                                  processing = "toil-rsem",
                                  unit="tpm",
                                  expect='datatable',
                                  read.from.redis = TRUE
                                  )
    expect_equal(expr0,NULL)
})
