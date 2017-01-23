library(futile.logger)
flog.appender(appender.file("test_setup.log"), name="log")

context("Setting up")

test_that("dbpath is set.", {
    expect_equal(dbpath[['gtex']], '/var/www/html/data/GTEx_V6.h5')
})

test_that("Redis server is running", {

    expect_error(rredis::redisConnect('blablahost'), "cannot open the connection", fixed=TRUE)
    expect_error(rredis::redisConnect('localhost'), NA)
    try ({
        rredis::redisSet("testvar", "this is a test variable")
        expect_equal(rredis::redisGet("testvar"), "this is a test variable")
    })

})

test_that("Load genelist into redis server", {
    fields = c("EnsemblID", "HGNC")
    start = proc.time()
    loadGeneData(cols=fields, write.to.redis = TRUE)
    runtime = proc.time() - start
    gDataRedis = rredis::redisGet('geneList')
    expect_equal(names(gDataRedis),fields)
    flog.info("loadGeneList to redis in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("redisGet('geneList'): ", head(gDataRedis), name="log", capture=TRUE)
})

test_that("Load limited genelist into redis server", {
    fields = c("EnsemblID", "HGNC")
    loadGeneData(cols=fields, write.to.redis = TRUE, genes=c('HS3ST1', 'HS3ST2', 'HS3ST3B1'))
    gDataRedis = rredis::redisGet('geneList')
    expect_equal(names(gDataRedis),fields)
    expect_equal(nrow(gDataRedis),3)
    flog.info("Load limited genelist into redis server", name="log")
    flog.info("redisGet('geneList'): ", head(gDataRedis), name="log", capture=TRUE)
})

test_that("Load gene list", {
    fields = c("EnsemblID", "HGNC")
    start = proc.time()
    loadGeneData(cols=fields, write.to.redis=FALSE)
    runtime = proc.time() - start
    gData = container$geneList
    expect_equal(names(gData), fields)
    flog.info("loadGeneList to R environment in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("container$geneList: ", gData, name="log", capture=TRUE)
})

test_that("Load gene sets into container", {
    loadGeneSets(write.to.redis = FALSE)
    getback = container$geneSets
    expect_equal(getback[[1]][['name']], 'Heparan sulfate 3-O-sulfotransferases')
})

test_that("Load gene sets into redis server", {
    loadGeneSets(write.to.redis = TRUE)
    getback = rredis::redisGet('geneSets')
    expect_equal(getback[[1]][['name']], 'Heparan sulfate 3-O-sulfotransferases')
})

test_that("Load sample metadata", {
    start = proc.time()
    loadSampleMetadata(write.to.redis=TRUE)
    runtime = proc.time() - start
    sampleData = rredis::redisGet('sampleMetadata')
    expect_equal(names(sampleData)[1], 'SAMPID')
    flog.info("loadSampleMetadata to Redis in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("redisGet('sampleMetadata'):", sampleData, name="log", capture=TRUE)

    # load again, only selected columns
    fields = c("SAMPID", "SMTS", "SMTSD")
    loadSampleMetadata(cols=fields, write.to.redis = TRUE)
    expect_equal(names(rredis::redisGet('sampleMetadata')), fields)

    # load all columns
    start = proc.time()
    loadSampleMetadata(write.to.redis=FALSE)
    runtime = proc.time() - start
    expect_equal(names(container$sampleMetadata)[1], 'SAMPID')
    flog.info("loadSampleMetadata in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("container$sampleMetadata:", container$sampleMetadata, name="log", capture=TRUE)

    # load again, only selected columns
    fields = c("SAMPID", "SMTS", "SMTSD")
    loadSampleMetadata(cols=fields)
    expect_equal(names(rredis::redisGet('sampleMetadata')), fields)
})



test_that("Load expression data", {

    ### define the subset
    samples = c("GTEX-UTHO-1226-SM-3GAEE",
                "GTEX-146FH-1726-SM-5QGQ2",
                "GTEX-QDT8-0126-SM-48TZ1",
                "GTEX-QCQG-1326-SM-48U24",
                "GTEX-WZTO-2926-SM-3NM9I",
                "GTEX-12WSB-0126-SM-59HJN",
                "GTEX-11VI4-0626-SM-5EQLO",
                "GTEX-T5JC-0526-SM-32PM7",
                "GTEX-RU1J-0426-SM-46MUK",
                "GTEX-1212Z-0226-SM-59HLF",
                "GTEX-11DZ1-2026-SM-5A5KG")
    genes.wVersion = c(
        "ENSG00000242268.2",
        "ENSG00000259041.1",
        "ENSG00000270112.3",
        "ENSG00000167578.16",
        "ENSG00000278814.1",
        "ENSG00000078237.5")
    genes = c(
        "ENSG00000242268",
        "ENSG00000259041",
        "ENSG00000270112",
        "ENSG00000167578",
        "ENSG00000278814",
        "ENSG00000078237")

    ### Load toil-rsem processed-data
    start = proc.time()
    loadExpressionData(genes, samples, db="gtex", processing="toil-rsem", unit="tpm", write.to.redis=FALSE,ctner=container)
    runtime = proc.time() - start
    m2 = get("/toil-rsem/gene/tpm/expressionMatrix", envir = container)
    expect_equal(genes, colnames(m2))
    expect_equal(samples, rownames(m2))

    start = proc.time()
    loadExpressionData(genes, samples, db="gtex", processing="toil-rsem", unit="tpm", write.to.redis=TRUE)
    runtime = proc.time() - start
    m2.redis = rredis::redisGet("/toil-rsem/gene/tpm/expressionMatrix")
    expect_equal(genes, colnames(m2.redis))
    expect_equal(samples, rownames(m2.redis))

    expect_equal(m2, m2.redis)

    flog.info("Load toil-rsem data subset in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("ls(container): ", ls(container), name="log", capture=TRUE)
    flog.info("expression matrix in container:", m2, name="log", capture=TRUE)
    flog.info("expression matrix in redis:", m2.redis, name="log", capture=TRUE)

    ### Load broad-processed data
    start = proc.time()
    loadExpressionData(genes, samples,db = "gtex", processing = "broad", unit = "fpkm", write.to.redis = FALSE, ctner = container)
    runtime = proc.time() - start
    m = get("/broad/gene/fpkm/expressionMatrix", envir = container)
    flog.info("Load broad data subset in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("ls(container): ", ls(container), name="log", capture=TRUE)
    flog.info("expression matrix:", m, name="log", capture=TRUE)

    # ### load everything, will take some time
    # start = proc.time()
    # loadExpressionData(db='gtex', processing = 'toil-rsem', unit = 'tpm',write.to.redis = TRUE)
    # runtime = proc.time() - start
    # m = rredis::redisGet('/toil-rsem/gene/tpm/expressionMatrix')
    # # expect_equal(ncol(m), 60498)
    # # expect_equal(nrow(m), 7863)
    # flog.info("loadExpressionMatrix in: %6.4f seconds.", runtime['elapsed'], name="log")
    # # flog.info("container$expressionMatrix:", str(container$expressionMatrix), name="log", capture=TRUE)
})
