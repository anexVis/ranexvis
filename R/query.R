# This will make available
# (list) dbpath
# (environment) container
data(sysdata,envir=environment())

#' getGeneList
#'
#' @param db the data set of interest: "gtex" is the only option available for now. For future development, db should be genomic feature annotation, such as genecode, or ensembl
#' @param cols a vector of column names to be retrieved: c('HGNC', 'EnsembleID', 'Description'), or NULL to retrieve all columns
#' @param expect a string specifying the output format: "json" for a json object (from jsonlite::toJSON), otherwise results in a data.table
#' @export
getGeneList <- function(db="gtex",cols=c('EnsemblID', 'HGNC'), expect='json',read.from.redis=TRUE) {
    if (read.from.redis) output = redisOpenGetClose("geneList")
    else output = get("geneList", envir = container)

    if (expect == 'json') {
        return(jsonlite::toJSON(output))
    } else {
        return(output)
    }
}

# getGeneSets <- function(db="gtex", cols=c("EnsemblID", "HGNC"), expec='json', read.from.redis=TRUE) {
#     if (read.from.redis) output = rredis::redisGet("geneSets")
#     else output = get("geneSets", envir=container)
#     if (expect == 'json') {
#         return(jsonlite::toJSON(output))
#     } else {
#         return(output)
#     }
# }

#' getSampleGroupingList
#'
#' @param db the data set of interest: "gtex" is the only option available for now. For future development, db should be genomic feature annotation, such as genecode, or ensembl
#' @param grouping available values are: 'SMTS' (tissue type), 'SMTSD' (sampled site), 'SMUBRID' (Uberon ID)
#' @export
getSampleGroupingList <- function(db="gtex", grouping="SMTS", expect='json',read.from.redis=TRUE) {
    if (read.from.redis) allmeta = redisOpenGetClose('sampleMetadata')
    else  allmeta = get("sampleMetadata", envir=container)
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
#' @param db database, default="gtex"
#' @param cols selected columns in metadata. Default value: c("SAMPID", "SMTS"). Available columns:
#' "SAMPID"     "SMATSSCR"   "SMNABTCH"   "SMNABTCHT"  "SMNABTCHD"  "SMGEBTCH"   "SMCENTER"   "SMPTHNTS"   "SMRIN"
#' "SMTS"       "SMTSD"      "SMUBRID"    "SMTSPAX"    "SMTSTPTREF" "SMAFRZE"
#' @param expect output format. Available values: 'json', 'datatable'
#' @export
getSampleMetadata <- function(db="gtex", cols=c("SAMPID", "SMTS"), expect="json", read.from.redis=TRUE) {
    if (read.from.redis)  allmeta = redisOpenGetClose("sampleMetadata")
    else allmeta = get("sampleMetadata", envir=container)
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
getGeneExpressionMatrix  <- function(genes, sampleGroups, sampleGrouping = "SMTS", db = "gtex", processing="toil-rsem", unit="tpm", expect="json",read.from.redis=TRUE) {
    sampleMeta = getSampleMetadata(db, cols=c("SAMPID", sampleGrouping), expect='datatable', read.from.redis=read.from.redis)

    path2dataset = paste("/", processing, "/gene/", unit, sep="")
    if (read.from.redis) fullExprMatrix =  redisOpenGetClose(paste(path2dataset, "expressionMatrix", sep="/"))
    else fullExprMatrix =  get(paste(path2dataset, "expressionMatrix", sep="/"), envir = container)

    genesInMatrix = colnames(fullExprMatrix)

    samplesInMatrix= data.table::data.table(rownames(fullExprMatrix))
    names(samplesInMatrix) = c("SAMPID")
    sampleList = merge(samplesInMatrix, sampleMeta, by = "SAMPID", all.x=TRUE,sort=FALSE)

    colidx = which(genesInMatrix %in% removeEnsemblVersion(genes))
    rowidx = which(sampleList[[sampleGrouping]] %in% sampleGroups)
    if (length(colidx) == 0 || length(rowidx) == 0) {
        message("No matching record found.")
        return(NULL)
    }

    tryCatch ({
        exprMatrix = fullExprMatrix[rowidx,colidx]
        if (expect=='datatable' || expect=='dt')  return(exprMatrix)
        else return(jsonlite::toJSON(exprMatrix))
    }, error = function(e) {
        print(e)
    })
}

#' Return ready-to-plot expression data for 2 genes
#' (This is a thin wrapper of getExpressionMatrix)
#'
getScatterData <- function(x,y, sampleGroups, sampleGrouping = "SMTS", db = "gtex", processing="toil-rsem", unit="tpm", read.from.redis=TRUE) {
    expr = getGeneExpressionMatrix(genes = c(x,y),
                                        sampleGroups=sampleGroups,
                                        sampleGrouping=sampleGrouping,
                                        db = db,
                                        processing=processing,
                                        unit = unit,
                                        expect='datatable',
                                        read.from.redis=FALSE
    )
    labels = ensembl2hgnc(removeEnsemblVersion(c(x,y)))
    names(expr)[1:2] = c('x', 'y')
    return(jsonlite::toJSON(list(xlabel=labels[1], ylabel=labels[2],data=expr)))
}
