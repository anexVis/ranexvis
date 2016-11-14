data(sysdata,envir=environment())


#' Load several data sets into the memory and store it in the globally available `container`
#'
setup <- function(genes=NULL, samples=NULL, write.to.redis=TRUE) {
    loadGeneData(write.to.redis=write.to.redis, ctner=container)
    loadSampleMetadata(swrite.to.redis=write.to.redis, ctner=container)
    loadExpressionData(genes=genes, samples=samples,write.to.redis=write.to.redis, ctner=container) # load everything will take about 30sec
    message("Finished loading all data to redis server")
}

#' Load the list of genes
#'
#' @param db default:'gtex', only option for now
#' @param cols default: c('EnsemblID', 'HGNC'), the fields to retrieve
#' @param write.to.redis default: TRUE. Write the data to redis if TRUE. Will invalidate `ctner` option.
#' @param ctner the environment to store the resulted gene list. Will be ignored if write.to.redis=TRUE
loadGeneData <- function(db="gtex",cols=c('EnsemblID', 'HGNC'), write.to.redis = TRUE, ctner=container) {
    path2dataset = paste("/genes", cols[1], sep="/")
    output = data.table::data.table(rhdf5::h5read(dbpath[[db]],path2dataset))
    for (column in cols[2:length(cols)]) {
        path2dataset = paste("/genes", column, sep="/")
        output[[column]] = as.character(rhdf5::h5read(dbpath[[db]], path2dataset))
    }
    names(output) = cols
    if (write.to.redis) {
        tryCatch(
            invisible(rredis::redisSet('geneList',output)),
            error = function(e){
                warning("Redis server not connected.")
            }
        )
    } else {
        invisible(assign("geneList", output, envir=ctner))
    }
}

loadSampleMetadata <- function(db='gtex', cols=NULL, write.to.redis = TRUE, ctner=container) {
    path2dataset = paste("/metadata/sample")
    rhdf5::H5close()    # does it affect h5 object handles of other process/user?
    output = data.table::data.table(rhdf5::h5read(dbpath[[db]], path2dataset))
    if (!is.null(cols)) {
        output = output[,cols]
    }
    if (write.to.redis) {
        tryCatch(
        invisible(rredis::redisSet('sampleMetadata',output)),
         error = function(e){
                warning("Redis server not connected.")
         })

    } else {
        invisible(assign("sampleMetadata", output, envir=ctner))
    }
}

loadExpressionData <- function(genes=NULL, samples=NULL,db = "gtex", processing="toil-rsem", unit="tpm", write.to.redis = TRUE, ctner=container) {
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
            if (write.to.redis) {
                invisible(checkConnectionAndSet(paste0(path2dataset, "/expressionMatrix"), NULL))
            } else {
                invisible(assign(paste0(path2dataset, "/expressionMatrix"), NULL,envir = ctner))
            }
        } else if (length(colidx) < length(genes)) {
            message(paste("Only", length(colidx), "in", length(genes), "are matched."))
        }

    }
    if (!is.null(samples)) {
        rowidx = which(sampleList %in% samples)
        if (length(rowidx) == 0) {
            message("No matching record found.")
            if (write.to.redis) {
                invisible(checkConnectionAndSet(paste0(path2dataset, "/expressionMatrix"), NULL))
            } else {
                invisible(assign(paste0(path2dataset,"/expressionMatrix"), NULL,envir = ctner))
            }
        }
    }

    # The subsetting of HDF5 is puzzling!
    # row-col seem to be switched between R and HDFView
    tryCatch ({
        exprMatrix = data.table::data.table(rhdf5::h5read(dbpath[[db]],
                                                            path2dataset,
                                                            index=list(rowidx,colidx)))
        if (!is.null(genes)) colnames(exprMatrix) = geneList[geneList %in% removeEnsemblVersion(genes)]
        else colnames(exprMatrix) = geneList

        if (!is.null(samples)) rownames(exprMatrix) = sampleList[sampleList %in% samples]
        else rownames(exprMatrix) = makeUniqueNames(sampleList)

        if (write.to.redis) {
            invisible(checkConnectionAndSet(paste0(path2dataset, "/expressionMatrix"),as.matrix( exprMatrix)))
        } else {
            invisible(assign(paste0(path2dataset, "/expressionMatrix"), as.matrix(exprMatrix),envir = ctner))
        }
    }, error = function(e) {
        print(e)
    })
}

checkConnectionAndSet <- function(key,value) {
    tryCatch(rredis::redisSet(key,value),
         error = function(e) {
             print(e)
             message("Value to set: ", str(value))
             stop("Value not set for key ", key)
    })
}
