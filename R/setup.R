data(sysdata,envir=environment())

container = new.env(parent= emptyenv())

#' Load several data sets into the memory
#'
setup <- function() {
    loadGeneData()
    loadSampleMetadata()
    loadExpressionData() # load everything will take about 30sec

}

loadGeneData <- function(db="gtex",cols=c('EnsemblID', 'HGNC')) {
    path2dataset = paste("/genes", cols[1], sep="/")
    output = data.table::data.table(rhdf5::h5read(dbpath[[db]],path2dataset))
    for (column in cols[2:length(cols)]) {
        path2dataset = paste("/genes", column, sep="/")
        output[[column]] = as.character(rhdf5::h5read(dbpath[[db]], path2dataset))
    }
    names(output) = cols
    invisible(assign("geneList", output, envir=container))
}

loadSampleMetadata <- function(db='gtex', cols=NULL) {
    path2dataset = paste("/metadata/sample")
    rhdf5::H5close()    # does it affect h5 object handles of other process/user?
    output = data.table::data.table(rhdf5::h5read(dbpath[[db]], path2dataset))
    if (!is.null(cols)) {
        output = output[,cols]
    }
    invisible(assign("sampleMetadata", output, envir=container))
}

loadExpressionData <- function(genes=NULL, samples=NULL,db = "gtex", processing="toil-rsem", unit="tpm") {
    path2dataset = paste("/", processing, "/gene/", unit, sep="")
    path2geneId   = paste("/", processing, "/gene/EnsemblID", sep="")
    path2sampleId = paste("/", processing, "/gene/SampleID", sep="")
    rhdf5::H5close()
    geneList = removeEnsemblVersion(readCharacterArray(dbpath[[db]], path2geneId)) # be careful about different version of gencode in each dataset
    sampleList = readCharacterArray(dbpath[[db]],path2sampleId)

    colidx = NULL
    rowidx = NULL
    if (!is.null(genes)) {
        colidx = which(geneList %in% removeEnsemblVersion(genes))
        if (length(colidx) == 0) {
            message("No matching record found.")
            invisible(assign(paste0(path2dataset, "expressionMatrix"), NULL,envir = container))
        }
    }
    if (!is.null(samples)) {
        rowidx = which(sampleList %in% samples)
        if (length(rowidx) == 0) {
            message("No matching record found.")
            invisible(assign(paste0(path2dataset,"/expressionMatrix"), NULL,envir = container))
        }
    }

    # The subsetting of HDF5 is puzzling!
    # row-col seem to be switched between R and HDFView
    tryCatch ({
        exprMatrix = data.table::data.table(rhdf5::h5read(dbpath[[db]],
                                                            path2dataset,
                                                            index=list(rowidx,colidx)))
        # colnames(exprMatrix) = genes
        # rownames(exprMatrix) = samples
        print(paste0(path2dataset, "/expressionMatrix"))
        invisible(assign(paste0(path2dataset, "/expressionMatrix"), exprMatrix,envir = container))
    }, error = function(e) {
        print(e)
    })
}
