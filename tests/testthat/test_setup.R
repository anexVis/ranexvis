library(futile.logger)
flog.appender(appender.file("test_setup.log"), name="log")

context("Setting up")
ctner = get("container")

test_that("Load gene list", {
    fields = c("EnsemblID", "HGNC")
    start = proc.time()
    loadGeneData(cols=fields)
    runtime = proc.time() - start
    # ctner = get("container")
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
    expect_equal(names(ctner$sampleMetadata)[1], 'SAMPID')
    flog.info("loadSampleMetadata in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("container$sampleMetadata:", ctner$sampleMetadata, name="log", capture=TRUE)

    # load again, only selected columns
    fields = c("SAMPID", "SMTS", "SMTSD")
    loadSampleMetadata(cols=fields)
    expect_equal(names(ctner$sampleMetadata), fields)
})

test_that("Load expression data", {


    ### load a subset
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
    genes = c(
        "ENSG00000242268.2",
        "ENSG00000259041.1",
        "ENSG00000270112.3",
        "ENSG00000167578.16",
        "ENSG00000278814.1",
        "ENSG00000078237.5")

    start = proc.time()
    loadExpressionData(genes, samples)
    runtime = proc.time() - start
    expect_equal(ncol(ctner$expressionMatrix), length(genes))
    expect_equal(nrow(ctner$expressionMatrix), length(samples))
    flog.info("loadExpressionMatrix subset in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("container$expressionMatrix:", ctner$expressionMatrix, name="log", capture=TRUE)


    ### load everything, will take some time
    start = proc.time()
    loadExpressionData()
    runtime = proc.time() - start
    expect_true(ncol(ctner$expressionMatrix) > 1)
    expect_true(nrow(ctner$expressionMatrix) > 1)
    flog.info("loadExpressionMatrix in: %6.4f seconds.", runtime['elapsed'], name="log")
    flog.info("container$expressionMatrix:", ctner$expressionMatrix, name="log", capture=TRUE)
})


