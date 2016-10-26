context("Analysis methods")

test_that("Co-expression calculation.", {
    selectedEnsembls = c("ENSG00000176022.3", "ENSG00000027847.9", "ENSG00000002587.5", "ENSG00000125430.4") # B3GALT6, B4GALT7
    expr = getGeneExpressionMatrix(genes = selectedEnsembls,
                                   sampleGroups=c("Brain"),
                                   sampleGrouping="SMTS",
                                   db = "gtex",
                                   processing = "toil-rsem",
                                   unit="tpm"
                                   )
    expect_equal(coexpression(expr), cor(expr,method="pearson")) # these two genes are known to be highly correlated
    expect_equal(coexpression(expr, method="spearman"), cor(expr,method="spearman")) # these two genes are known to be highly correlated
    # expect_equal(cor(expr,method="spearman")[1,1], 1)
    # expect_true(cor(expr,method="sp")[1,2] > 0.5)
})
