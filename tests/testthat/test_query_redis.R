library(futile.logger)
context("Querying from redis server")
flog.appender(appender.file("test_query_redis.log"), name="log")

# connecting to redis was done on startup (should it be done with individual queries instead?)
# setup() was called upon startup

test_that("Redis server is online and connected.", {
    expect_error(rredis::redisConnect('localhost'), NA)
    pilot_genes = c('XYLT1', 'XYLT2','B4GALT7', 'B3GALT6',
                    'B3GAT3', 'EXTL2', 'EXTL3', 'EXTL1',
                    'EXT1', 'EXT2', 'NDST1', 'NDST2',
                    'NDST3', 'NDST4', 'GLCE', 'HS2ST1',
                    'HS6ST1', 'HS6ST2', 'HS6ST3', 'HS3ST1',
                    'HS3ST2', 'HS3ST3A1', 'HS3ST3B1', 'HS3ST4',
                    'HS3ST5', 'HS3ST6')
    expect_message(setup(genes=hgnc2ensembl(pilot_genes), write.to.redis = TRUE),
                   "Finished loading data to redis server")
    tmp = rredis::redisGet('/toil-rsem/gene/tpm/expressionMatrix')
    expect_false(is.null(tmp))
    flog.info("expressionMatrix in redis server: ", dim(tmp), name='log', capture = TRUE)
    expect_equal(length(pilot_genes), ncol(tmp))
    expect_true(nrow(tmp) > 500)
})


test_that("Retrieval of gene list.", {
    # expect_warning(genelist = getGeneList(expect='df'),'HDF5 on bit64conversion')
    # expect_warning(jsonlist = getGeneList(db='gtex', expect='json'), 'HDF5 on bit64conversion')
    genelist = getGeneList("gtex", c("EnsemblID", "HGNC"), expect="df")
    expect_equal(names(genelist) , c("EnsemblID", "HGNC"))
})

test_that("Retrieval of sample grouping list.", {
   sampleGroupingList = getSampleGroupingList(expect='df')
})

test_that("Retrieval of sample metadata.", {
   sampleMetadata = getSampleMetadata("gtex", cols=c("SAMPID", "SMTS"), expect="dt")
   expect_equal(ncol(sampleMetadata), 2)
})

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
                                  expect='datatable'
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
                                  expect='datatable'
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
                                  expect='datatable'
                                  ), "No matching record found.")
})

