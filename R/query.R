#' getGeneList
#'
#' @param db the large experiment of interest: "gtex" is the only option available for now.
#' @param workflow the processing workflow to obtain normalized RNA expression level, available values are "broad", "toil-rsem", "toil-kallisto"
#' @param unit the unit for measuring abundance: `tpm` or `fpkm`
#' @param cols a vector of column names to be retrieved: c('HGNC', 'EnsembleID', 'Description'), or NULL to retrieve all columns
#' @param expect a string specifying the output format: "json" for a json object (from jsonlite::toJSON), otherwise results in a data.frame
#' @export
#' @import rhdf5
getGeneList <- function(db="gtex",cols=c('EnsemblID', 'HGNC'), expect='json') {
    get("dbpath")
    output = list()
    for (column in cols) {
        output[[column]] = h5read(dbpath[[db]],name = paste("/genes",column, sep="/"))
    }

    if (expect == 'json') {
        return(toJSON(data.frame(output)))
    } else {
        return(data.frame(output))
    }
}

