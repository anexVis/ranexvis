#' Utility function to read a 1D array correctly from HDF5 file
#'
#' @param file path to HDF5 file
#' @param path path to the dataset within the HDF5
#' @param colname if NULL, this will return a character vector. Otherwise return a data.frame with 1 column named `colname`.
readCharacterArray <- function(file, path, colname=NULL) {
    charray = as.character(rhdf5::h5read(file, path))
    if (is.null(colname)) {
        return(charray)
    } else {
        output = list()
        output[[colname]] = charray
        return(data.frame(output))
    }
}

removeEnsemblVersion <- function(x) {
    return(gsub("(ENSG\\d+)\\.(\\d+)", "\\1", x))
}
