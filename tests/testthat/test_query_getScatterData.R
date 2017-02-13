library(futile.logger)
context("Molding scatter plot data from redis server")
flog.appender(appender.file("test_query_getScatterData.log"), name="log")

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

test_that("Retrieval of scatter data with sample meta data", {
    expr_json = getScatterData(x = "ENSG00000176022.3",
                               y = "ENSG00000027847.9",
                               sampleGroups="Bladder",
                               sampleGrouping="SMTSD",
                               sampleMetaFields = c("SMTS", "SMTSD", "GENDER", "RACE"),
                               db = "gtex",
                               processing = "toil-rsem",
                               unit="tpm",
                               read.from.redis=TRUE
    )
   flog.info("Scatterplot data with sample meta data, in JSON:", jsonlite::prettify(expr_json), name='log', capture=TRUE)
})

test_that("Retrieval of scatter data without sample meta data", {
    expr_json = getScatterData(x = "ENSG00000176022.3",
                               y = "ENSG00000027847.9",
                               sampleGroups="Bladder",
                               sampleGrouping="SMTSD",
                               db = "gtex",
                               processing = "toil-rsem",
                               unit="tpm",
                               read.from.redis=TRUE
    )
    expr_obj1 = jsonlite::fromJSON(expr_json)

    # Switch the order
    expr_json = getScatterData(y = "ENSG00000176022.3",
                               x = "ENSG00000027847.9",
                               sampleGroups="Bladder",
                               sampleGrouping="SMTSD",
                               db = "gtex",
                               processing = "toil-rsem",
                               unit="tpm",
                               read.from.redis=TRUE
    )

    expr_obj2 = jsonlite::fromJSON(expr_json)
    expect_equal(expr_obj1$dimensions[['name']], c('x', 'y'))
    expect_equal(expr_obj1$dimensions[['label']], c('B3GALT6', 'B4GALT7'))

    expect_equal(expr_obj2$dimensions[['name']], c('x', 'y'))
    expect_equal(expr_obj2$dimensions[['label']], c('B4GALT7', 'B3GALT6'))

    expect_equal(expr_obj1$data[['x']], expr_obj2$data[['y']])
    expect_equal(expr_obj1$data[['y']], expr_obj2$data[['x']])

    flog.info("Scatterplot data without sample meta data, in JSON:", jsonlite::prettify(expr_json), name='log', capture=TRUE)
})

test_that("Retrieval of scatter data of identical pair", {
    expr_json = getScatterData(x = "ENSG00000176022.3",
                               y = "ENSG00000176022.3",
                               sampleGroups="Bladder",
                               sampleGrouping="SMTSD",
                               db = "gtex",
                               processing = "toil-rsem",
                               unit="tpm",
                               read.from.redis=TRUE
    )
    expr_obj = jsonlite::fromJSON(expr_json)
    expect_equal(expr_obj$dimensions[['name']], 'x')
    expect_equal(expr_obj$dimensions[['label']], 'B3GALT6')
    flog.info("Scatterplot data of identical pair, in JSON:", jsonlite::prettify(expr_json), name='log', capture=TRUE)
})

test_that("Retrieval of scatter data of identical pair with metadata", {
    expr_json = getScatterData(x = "ENSG00000176022",
                               y = "ENSG00000176022",
                               sampleGroups="Bladder",
                               sampleGrouping="SMTSD",
                               sampleMetaFields = c("SMTSD"),
                               db = "gtex",
                               processing = "toil-rsem",
                               unit="tpm",
                               read.from.redis=TRUE
    )
    expr_obj = jsonlite::fromJSON(expr_json)
    expect_equal(expr_obj$dimensions[['name']], c('x', 'SMTSD'))
    flog.info("Scatterplot data of identical pair with metadata, in JSON:", jsonlite::prettify(expr_json), name='log', capture=TRUE)
})
