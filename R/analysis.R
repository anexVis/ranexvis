#' Evaluate pairwise co-expression using one of the methods: Pearson, Spearman, Mutual Information (to be implemented)
#'
#' @param x expression matrix (samples x genes)
#' @param method valid options: "pearson" (default), "spearman", "mi" (not yet implemented)
#'
#' For additional parameters such as `na.rm`, `use` refer to the documentation on `cor {stats}`.
coexpression <- function(x, method="pearson",...) {
    corrMatrix =
    switch( method,
            "pearson" = cor(x, method="pearson"),
            "spearman" = cor(x, method="spearman"),
            "mi" = NULL,
            cor(x, method="pearson"))
    return(corrMatrix)
}

#' Return ready-to-use json data for plotting of heatmap, including a matrix, a row dendrogram, and a column dendrogram
#' This is a convenient wrapping of rglyvis::coexpression and rvislib::heatmap.adjacency
#'
#' @export
coexpression.heatmap <- function(genes, sampleGroups, sampleGrouping, db,processing, unit,...) {
    exprMatrix = getGeneExpressionMatrix(genes, sampleGroups, sampleGrouping, db, processing, unit)
    colnames(exprMatrix) = genes
    corrMatrix = coexpression(exprMatrix, ...)
    return(rvislib::heatmap.adjacency(corrMatrix))
}

