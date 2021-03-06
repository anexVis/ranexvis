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
getGeneList <- function(db="gtex",cols=c('EnsemblID', 'HGNC'), withEnsemblVersion = TRUE, expect='json',read.from.redis=TRUE) {
    if (read.from.redis) output = redisOpenGetClose("geneList")
    else output = get("geneList", envir = container)

    if (!withEnsemblVersion && 'EnsemblID' %in% cols) {
        output[['EnsemblID']] = removeEnsemblVersion(output[['EnsemblID']])
    }

    if (expect == 'json') {
        return(jsonlite::toJSON(output))
    } else {
        return(output)
    }
}

#' getSampleGroupingList
#'
#' @param db the data set of interest: "gtex" is the only option available for now. For future development, db should be genomic feature annotation, such as genecode, or ensembl
#' @param grouping available values are: 'SMTS' (tissue type), 'SMTSD' (sampled site), 'SMUBRID' (Uberon ID). If not specified, 'SMTS' will be returned
#' @export
getSampleGroupingList <- function(db="gtex", grouping="SMTS", expect='json',read.from.redis=TRUE) {
    if (read.from.redis) allmeta = redisOpenGetClose('sampleMetadataPhenotype')
    else  allmeta = get("sampleMetadataPhenotype", envir=container)
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
#' @param cols selected columns in metadata. Default value: c("SAMPID", "SMTS"). Available columns includes both sample attributes such as:
#' "SAMPID"     "SMATSSCR"   "SMNABTCH"   "SMNABTCHT"  "SMNABTCHD"  "SMGEBTCH"   "SMCENTER"   "SMPTHNTS"   "SMRIN"
#' "SMTS"       "SMTSD"      "SMUBRID"    "SMTSPAX"    "SMTSTPTREF" "SMAFRZE"
#' and subject phenotypes such as
#' "AGE", "GENDER", "RACE", "ETHNCTY", "HGHT", "WGHT", "BMI"
#' @param expect output format. Available values: 'json', 'datatable'
#' @export
getSampleMetadata <- function(db="gtex", sampleIds=NULL , cols=c("SAMPID", "SMTS"), expect="json", read.from.redis=TRUE) {
    if (read.from.redis)  allmeta = redisOpenGetClose("sampleMetadataPhenotype")
    else allmeta = get("sampleMetadataPhenotype", envir=container)

    if (is.null(sampleIds))  output = allmeta[cols]
    else output = allmeta[allmeta$SAMPID %in% sampleIds,cols]

    returnData = switch (expect,
                         'json' = jsonlite::toJSON(output),
                         'datatable' = output,
                         output
    )
    return(returnData)
}



#' getSampleMetadataByGroup
#'
#' A shortcut for querying sample meta data by some grouping
#' @param sampleGroups String or string vector, name of the groups to be selected, e.g. c('Brain', 'Liver')
#' @param sampleGrouping [="SMTS"] available options: "SMTS", "SMTSD"
#' @param db [="gtex"]
#' @param cols [=c("SAMPID")]
#' @param expect [="json"]
#' @param read.from.redis [=TRUE]
#' @export
getSampleMetadataByGroup <- function(sampleGroups, sampleGrouping = "SMTS", db = "gtex", cols=c('SAMPID'), expect="json",read.from.redis=TRUE) {
    cols = union(cols,c("SAMPID", sampleGrouping))
    if (read.from.redis)  allmeta = redisOpenGetClose("sampleMetadataPhenotype")
    else allmeta = get("sampleMetadata", envir=container)

    output = allmeta[ allmeta[[sampleGrouping]] %in% sampleGroups,cols]

    returnData = switch (expect,
           'json' = jsonlite::toJSON(output),
           'datatable' = output,
           output
            )
    return(returnData)
}

#' Return a samples x genes expression matrix for correlation calculation
#'
#'@export
getGeneExpressionMatrix  <- function(genes, sampleGroups, sampleGrouping = "SMTS", sampleMetaFields=NULL,db = "gtex", processing="toil-rsem", unit="tpm", expect="json",read.from.redis=TRUE) {

    sampleMeta = getSampleMetadata(db, cols=union(c("SAMPID", sampleGrouping), sampleMetaFields), expect='datatable', read.from.redis=read.from.redis)

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
        exprDTable = data.table::data.table(exprMatrix)
        # 1-column matrix needs special treatment
        if (class(exprMatrix)[1] == 'numeric') {
            names(exprDTable)[1] = removeEnsemblVersion(genes)[1]
        }
        if (!is.null(sampleMetaFields)) {
            selectedSamples = (data.table::as.data.table(sampleList[rowidx, sampleMetaFields]))
            names(selectedSamples) = sampleMetaFields
            exprDTable= cbind.fast(exprDTable,selectedSamples)
        }
        if (expect=='datatable' || expect=='dt')  {
            return(exprDTable)
        } else return(jsonlite::toJSON(exprDTable, digits=I(8),force = TRUE))
    }, error = function(e) {
        stop(paste("Error in getGeneExpressionMatrix:", e))
    })
}

#' Return ready-to-plot expression data for 2 genes,
#' plus the metadata of each point, if specified
#' (This is a thin wrapper of getGeneExpressionMatrix)
#'
#' @export
getScatterData <- function(x,y, sampleGroups, sampleGrouping = "SMTS", sampleMetaFields=NULL,db = "gtex", processing="toil-rsem", unit="tpm", read.from.redis=TRUE) {
    if (x == y) genes = x
    else genes = c(x,y)
    expr = getGeneExpressionMatrix(genes = genes,
                                        sampleGroups=sampleGroups,
                                        sampleGrouping=sampleGrouping,
                                        sampleMetaFields=sampleMetaFields,
                                        db = db,
                                        processing=processing,
                                        unit = unit,
                                        expect='datatable',
                                        read.from.redis=read.from.redis
    )
    if (x == y) {
        names(expr)[1] = 'x'
        labels = c(ensembl2hgnc(removeEnsemblVersion(x)),sampleMetaFields)
        dimensions = data.frame(list(name=names(expr),label=labels))
    } else {
        ensemblIds = removeEnsemblVersion(c(x,y))
        labels = c(ensembl2hgnc(ensemblIds),sampleMetaFields)
        names(expr)[1:2] = c('x', 'y')[match(names(expr[1:2]), ensemblIds)]
        if (is.null(sampleMetaFields))
            dimensions = data.frame(list(name=c('x', 'y'),label=labels))
        else
            dimensions = data.frame(list(name=c('x', 'y', names(expr)[3:ncol(expr)]),label=labels))
    }
    return(jsonlite::toJSON(list(dimensions=dimensions,data=expr), digits=I(8)))
}


#' Return a array of gene set, each item having the format as in the example:
#' {name: 'HS synthetic enzymes', value: ['ENSG1', 'ENSG2',...]}
#'
#' @export
getGeneSets <- function(db="gtex", processing="toil-rsem", expect="json",read.from.redis=TRUE) {
    if (read.from.redis)
        output = redisOpenGetClose('geneSets')
    else
        output = geneSets

    if (expect=='json')
        return(jsonlite::toJSON(output,auto_unbox = TRUE))
    else
        return(output)
}
