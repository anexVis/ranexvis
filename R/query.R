data(sysdata,envir=environment())
ctner = get("container")

#' getGeneList
#'
#' @param db the data set of interest: "gtex" is the only option available for now. For future development, db should be genomic feature annotation, such as genecode, or ensembl
#' @param cols a vector of column names to be retrieved: c('HGNC', 'EnsembleID', 'Description'), or NULL to retrieve all columns
#' @param expect a string specifying the output format: "json" for a json object (from jsonlite::toJSON), otherwise results in a data.table
#' @export
getGeneList <- function(db="gtex",cols=c('EnsemblID', 'HGNC'), expect='json') {
    output = get("geneList", envir = container)
    if (expect == 'json') {
        return(jsonlite::toJSON(output))
    } else {
        return(output)
    }
}

#' getSampleGroupingList
#'
#' @param db the data set of interest: "gtex" is the only option available for now. For future development, db should be genomic feature annotation, such as genecode, or ensembl
#' @param grouping available values are: 'SMTS' (tissue type), 'SMTSD' (sampled site), 'SMUBRID' (Uberon ID)
#' @export
getSampleGroupingList <- function(db="gtex", grouping="SMTS", expect='json') {
    # TODO not sure how to read only one column. That would be more efficient
    path2dataset = paste("/metadata/sample")
    rhdf5::H5close()    # does it affect h5 object handles of other process/user?
    allmeta = data.table::data.table(rhdf5::h5read(dbpath[[db]], path2dataset))
    output = list()
    output[[grouping]] = unique(allmeta[[grouping]])
    if (expect == 'json') {
        return(jsonlite::toJSON(data.frame(output)))
    } else {
        return(data.frame(output))
    }
}

#' getSampleMetadata
#'
#' @param db
#' @param cols selected columns in metadata. Default value: c("SAMPID", "SMTS"). Available columns:
#' "SAMPID"     "SMATSSCR"   "SMNABTCH"   "SMNABTCHT"  "SMNABTCHD"  "SMGEBTCH"   "SMCENTER"   "SMPTHNTS"   "SMRIN"
#' "SMTS"       "SMTSD"      "SMUBRID"    "SMTSPAX"    "SMTSTPTREF" "SMAFRZE"
#' @param expect output format. Available values: 'json', 'datatable'
#' @export
getSampleMetadata <- function(db="gtex", cols=c("SAMPID", "SMTS"), expect="json") {
    path2dataset = paste("/metadata/sample")
    rhdf5::H5close()
    allmeta = data.table::data.table(rhdf5::h5read(dbpath[[db]], path2dataset))
    output = allmeta[cols]
    returnData = switch (expect,
           'json' = jsonlite::toJSON(output),
           'datatable' = output,
           output
            )
    return(returnData)
}

#' Return a samples x genes expression matrix for correlation calculation
#'
getGeneExpressionMatrix  <- function(genes, sampleGroups, sampleGrouping = "SMTS", db = "gtex", processing="toil-rsem", unit="tpm") {
    sampleMeta = getSampleMetadata(db, cols=c("SAMPID", sampleGrouping), expect="datatable")
    path2dataset = paste("/", processing, "/gene/", unit, sep="")
    path2geneId   = paste("/", processing, "/gene/EnsemblID", sep="")
    path2sampleId = paste("/", processing, "/gene/SampleID", sep="")
    rhdf5::H5close()
    geneList = removeEnsemblVersion(readCharacterArray(dbpath[[db]], path2geneId)) # be careful about different version of gencode in each dataset
    sampleList = merge(readCharacterArray(dbpath[[db]],path2sampleId, colname = "SAMPID"), sampleMeta, by = "SAMPID", all.x=TRUE)
    colidx = which(geneList %in% removeEnsemblVersion(genes))
    rowidx = which(sampleList[[sampleGrouping]] %in% sampleGroups)
    if (length(colidx) == 0 || length(rowidx) == 0) {
        message("No matching record found.")
        return(NULL)
    }

    # The subsetting of HDF5 is puzzling!
    # row-col seem to be switched between R and HDFView
    tryCatch ({
        exprMatrix = data.table::data.table(rhdf5::h5read(dbpath[[db]],
                                                            path2dataset,
                                                            index=list(rowidx,colidx)))
        colnames(exprMatrix) = genes
        return(exprMatrix)
    }, error = function(e) {
        print(e)
    })
}
