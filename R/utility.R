# TODO: integrate these data into local database
    # ensembl servers usually broken
    # HOTFIX: change host between main server and archive server
    # www.ensembl.org
    # oct2016.archive.ensembl.org
ensemblHost = 'oct2016.archive.ensembl.org'

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

makeUniqueNames <- function(x) {
    count = table(x)
    renamed = x
    j = 0
    for (i in 1:length(x)) {
        xi = x[i]
        if (count[xi] == 1) next
        else {
            if (j == count[xi]) j = 0
            j = j+1
            renamed[i] = (paste(xi, j, sep="."))
        }
    }
    return(renamed)
}

ensembl2entrez <- function(ensemblIDs) {
    request =paste('http://rest.ensembl.org/xrefs/id/',
                   ensemblID,
                   '?content-type=application/json;external_db=EntrezGene', sep='')
    tmp <- jsonlite::fromJSON(request)
    return(tmp$primary_id)
}

# hgnc2ensembl <- function(hgnc_symbol) {
#     # the var name will be made col name in data.table. keep it the same to merge
#     ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
#                                dataset = 'hsapiens_gene_ensembl',
#                                host='www.ensembl.org')
#     rep_ensembl = biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
#                         filter='hgnc_symbol', values=hgnc_symbol, mart=ensembl)
#     ids = data.table::data.table(hgnc_symbol)
#     names(ids) = c('hgnc_symbol')
#     tmp = merge(ids, rep_ensembl, by='hgnc_symbol', sort=FALSE)
#     return(tmp[,2])
# }

hgnc2ensembl <- function(hgnc_symbol) {
    # the var name will be made col name in data.table. keep it the same to merge

    ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset = 'hsapiens_gene_ensembl',
                               host=ensemblHost)
    rep_ensembl = biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'ucsc'),
                                 filter='hgnc_symbol', values=hgnc_symbol, mart=ensembl)
    ids = data.table::data.table(hgnc_symbol)
    names(ids) = c('hgnc_symbol')
    tmp = merge(ids, rep_ensembl, by='hgnc_symbol', sort=FALSE)
    # some genes have duplicate Ensembl ID, mapping to different regions in different assembly
    # further filtering on 'ucsc' to make sure that only the ID in the primary assembly is used
    tmp = subset(tmp, ucsc != "")
    return(unique(tmp[['ensembl_gene_id']]))
}


ensembl2hgnc.remote <- function(ensembl_id) {
    # the var name will be made col name in data.table. keep it the same to merge
    ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset = 'hsapiens_gene_ensembl',
                               host=ensemblHost)
    rep_ensembl = biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                        filter='ensembl_gene_id', values=ensembl_id, mart=ensembl)
    ids = data.table::data.table(ensembl_id)
    names(ids) = c('ensembl_gene_id')
    tmp = merge(ids, rep_ensembl, by='ensembl_gene_id', sort=FALSE)
    return(tmp[,2])
}

ensembl2hgnc <- function(ensembl_id) {
    glist = redisOpenGetClose('geneList')
    return(glist[removeEnsemblVersion(EnsemblID)==removeEnsemblVersion(ensembl_id),HGNC])
}


#' Wrapping around redisGet to use with OpenCPU
#' Since each function call via OpenCPU does not
#' have access to redis connection created by other calls,
#' this function is used to create a short-lived connection
#' that will be closed immediately after value is retrieved
#'
redisOpenGetClose <- function(key) {
    rredis::redisConnect()
    val = rredis::redisGet(key)
    rredis::redisClose()
    return(val)
}

redisOpenSetClose <- function(key, value) {
    rredis::redisConnect()
    rredis::redisSet(key,value)
    rredis::redisClose()
}


#' fast cbind for data.table (about 30 times faster than cbind or cbind2)
#'
cbind.fast <- function(...) {
    x = c(...)
    data.table::setattr(x, "class", c("data.table", "data.frame"))
    ans = .Call(data.table:::Calloccolwrapper, x, max(100L, ncol(x) + 64L), FALSE)
    .Call(data.table:::Csetnamed, ans, 0L)
}
