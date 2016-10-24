data(sysdata,envir=environment())


#' getGeneList
#'
#' @param db the data set of interest: "gtex" is the only option available for now. For future development, db should be genomic feature annotation, such as genecode, or ensembl
#' @param cols a vector of column names to be retrieved: c('HGNC', 'EnsembleID', 'Description'), or NULL to retrieve all columns
#' @param expect a string specifying the output format: "json" for a json object (from jsonlite::toJSON), otherwise results in a data.table
#' @export
getGeneList <- function(db="gtex",cols=c('EnsemblID', 'HGNC'), expect='json') {
    path2dataset = paste("/genes", cols[1], sep="/")
    output = tryCatch({
        data.table::data.table(rhdf5::h5read(dbpath[[db]],path2dataset))
    }, error = function(e) {
        print(paste("Trying to read dataset", path2dataset, "from h5 file",dbpath[[db]]))
        print(e)
    })
    for (column in cols[2:length(cols)]) {
        path2dataset = paste("/genes", column, sep="/")
        output[[column]] = as.character(rhdf5::h5read(dbpath[[db]], path2dataset))
    }
    names(output) = cols
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
    allmeta = data.table::data.table(rhdf5::h5read(dbpath[[db]], path2dataset))
    output = list()
    output[[grouping]] = unique(allmeta[[grouping]])
    if (expect == 'json') {
        return(jsonlite::toJSON(data.frame(output)))
    } else {
        return(data.frame(output))
    }
}
