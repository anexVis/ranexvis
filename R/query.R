#' getGeneList
#'
#' @param db the large experiment of interest: "gtex" is the only option available for now. For future development, db should be genomic feature annotation, such as genecode, or ensembl
#' @param cols a vector of column names to be retrieved: c('HGNC', 'EnsembleID', 'Description'), or NULL to retrieve all columns
#' @param expect a string specifying the output format: "json" for a json object (from jsonlite::toJSON), otherwise results in a data.frame
#' @export
#' @import rhdf5
#' @import jsonlite
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
