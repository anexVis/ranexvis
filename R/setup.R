data(sysdata,envir=environment())


#' Load several data sets into the memory and store it in the globally available `container`
#'
setup <- function() {
    loadGeneData(ctner=container)
    loadSampleMetadata()
    loadExpressionData() # load everything will take about 30sec

}

loadGeneData <- function(db="gtex",cols=c('EnsemblID', 'HGNC'), ctner=container) {
    path2dataset = paste("/genes", cols[1], sep="/")
    output = data.table::data.table(rhdf5::h5read(dbpath[[db]],path2dataset))
    for (column in cols[2:length(cols)]) {
        path2dataset = paste("/genes", column, sep="/")
        output[[column]] = as.character(rhdf5::h5read(dbpath[[db]], path2dataset))
    }
    names(output) = cols
    invisible(assign("geneList", output, envir=ctner))
}

loadSampleMetadata <- function(db='gtex', cols=NULL, ctner=container) {
    path2dataset = paste("/metadata/sample")
    rhdf5::H5close()    # does it affect h5 object handles of other process/user?
    output = data.table::data.table(rhdf5::h5read(dbpath[[db]], path2dataset))
    if (!is.null(cols)) {
        output = output[,cols]
    }
    invisible(assign("sampleMetadata", output, envir=ctner))
}

loadExpressionData <- function(genes=NULL, samples=NULL,db = "gtex", processing="toil-rsem", unit="tpm", ctner=container) {
    path2dataset = paste("/", processing, "/gene/", unit, sep="")
    path2geneId   = paste("/", processing, "/gene/EnsemblID", sep="")
    path2sampleId = paste("/", processing, "/gene/SampleID", sep="")
    rhdf5::H5close()
    # be careful about different version of gencode in each dataset
    # the query usually does not care about specific version of an Ensembl ID
    geneList = removeEnsemblVersion(readCharacterArray(dbpath[[db]], path2geneId))
    sampleList = readCharacterArray(dbpath[[db]],path2sampleId)

    colidx = NULL
    rowidx = NULL
    if (!is.null(genes)) {
        colidx = which(geneList %in% removeEnsemblVersion(genes))
        if (length(colidx) == 0) {
            message("No matching record found.")
            invisible(assign(paste0(path2dataset, "expressionMatrix"), NULL,envir = ctner))
        } else if (length(colidx) < length(genes)) {
            message(paste("Only", length(colidx), "in", length(genes), "are matched."))
        }

    }
    if (!is.null(samples)) {
        rowidx = which(sampleList %in% samples)
        if (length(rowidx) == 0) {
            message("No matching record found.")
            invisible(assign(paste0(path2dataset,"/expressionMatrix"), NULL,envir = ctner))
        }
    }
c
    # The subsetting of HDF5 is puzzling!
    # row-col seem to be switched between R and HDFView
    tryCatch ({
        exprMatrix = data.table::data.table(rhdf5::h5read(dbpath[[db]],
                                                            path2dataset,
                                                            index=list(rowidx,colidx)))
        colnames(exprMatrix) = geneList[geneList %in% removeEnsemblVersion(genes)]
        rownames(exprMatrix) = sampleList[sampleList %in% samples]
        invisible(assign(paste0(path2dataset, "/expressionMatrix"), exprMatrix,envir = ctner))
    }, error = function(e) {
        print(e)
    })
}
