library(stringr)
library(futile.logger)

# setup()

flog.appender(appender.file("test_analysis.log"), 'log')
context("Analysis methods")
# selectedEnsembls = c("ENSG00000176022.3", "ENSG00000027847.9", "ENSG00000002587.5", "ENSG00000125430.4")
selectedEnsembls = c("ENSG00000176022", "ENSG00000027847", "ENSG00000002587", "ENSG00000125430")
ensembl2genes = c('B3GALT6', 'B4GALT7', 'HS3ST1', 'HS3ST3B1')
expr = getGeneExpressionMatrix(genes = selectedEnsembls,
                                   sampleGroups=c("Brain"),
                                   sampleGrouping="SMTS",
                                   db = "gtex",
                                   processing = "toil-rsem",
                                   unit="tpm",
                                   expect = 'dt'
                                   )
test_that("Co-expression calculation.", {
    flog.info("Co-expression calculation", name="log")
    pcorr = coexpression(expr)
    scorr = coexpression(expr, method="spearman")
    expect_equal(pcorr, cor(expr,method="pearson")) # these two genes are known to be highly correlated
    expect_equal(scorr, cor(expr,method="spearman")) # these two genes are known to be highly correlated
    flog.info("Pearson correlation", pcorr, name="log", capture=TRUE)
    flog.info("Spearman correlation", scorr, name="log", capture=TRUE)
})

test_that("Heatmap data generator", {
    hmdata = coexpression.heatmap(selectedEnsembls, "Brain",sampleGrouping = "SMTS", db = "gtex", processing="toil-rsem", unit="tpm")
    expect_equal(length(hmdata), 2)
    expect_equal(class(hmdata[1]), "character")

    flog.info("Matrix data ", jsonlite::prettify(hmdata[1]) , name="log", capture=TRUE)
    matrixData = jsonlite::fromJSON(hmdata[1])

})
