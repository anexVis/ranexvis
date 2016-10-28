library(stringr)
library(futile.logger)

flog.appender(appender.file("test_analysis.log"), 'log')
context("Analysis methods")
selectedEnsembls = c("ENSG00000176022.3", "ENSG00000027847.9", "ENSG00000002587.5", "ENSG00000125430.4") # B3GALT6, B4GALT7
    expr = getGeneExpressionMatrix(genes = selectedEnsembls,
                                   sampleGroups=c("Brain"),
                                   sampleGrouping="SMTS",
                                   db = "gtex",
                                   processing = "toil-rsem",
                                   unit="tpm"
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
    expect_equal(class(hmdata[0]), "character")
    flog.info("Heatmap data returned", hmdata, name="log", capture=TRUE)
})
