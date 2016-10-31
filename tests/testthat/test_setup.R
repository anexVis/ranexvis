library(futile.logger)
flog.appender(appender.file("test_setup.log"), name="log")

context("Setting up")

test_that("Load gene list", {
    fields = c("EnsemblID", "HGNC")
    start = proc.time()
    loadGeneData(cols=fields)
    runtime = proc.time() - start
    ctner = get("container")
    gData = ctner$geneList
    expect_equal(names(gData), fields)
    flog.info("loadGeneList in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("container$geneList: ", gData, name="log", capture=TRUE)
})

test_that("Load sample metadata", {
    # load all columns
    start = proc.time()
    loadSampleMetadata()
    runtime = proc.time() - start
    ctner = get("container")
    expect_equal(names(ctner$sampleMetadata)[1], 'SAMPID')
    flog.info("loadSampleMetadata in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("container$sampleMetadata:", ctner$sampleMetadata, name="log", capture=TRUE)

    # load again, only selected columns
    fields = c("SAMPID", "SMTS", "SMTSD")
    loadSampleMetadata(cols=fields)
    expect_equal(names(ctner$sampleMetadata), fields)
})
