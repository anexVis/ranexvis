pkg.env <- new.env(parent=emptyenv())
.onLoad <- function(libname, pkgname) {
    # open redis connection
    rredis::redisConnect('localhost',returnRef = TRUE)

    ### for development only
    # loading full set of expression data is not possible with the current implementation and setup:
    # redis server has a limit of 2GB-value size, and redis cluster has to be installed on multiple physical node
    loadGeneSets(write.to.redis = TRUE)
    # pilot_genes = c('XYLT1', 'XYLT2','B4GALT7', 'B3GALT6',
    #                 'B3GAT3', 'EXTL2', 'EXTL3', 'EXTL1',
    #                 'EXT1', 'EXT2', 'NDST1', 'NDST2',
    #                 'NDST3', 'NDST4', 'GLCE', 'HS2ST1',
    #                 'HS6ST1', 'HS6ST2', 'HS6ST3', 'HS3ST1',
    #                 'HS3ST2', 'HS3ST3A1', 'HS3ST3B1', 'HS3ST4',
    #                 'HS3ST5', 'HS3ST6')

    pilot_genes = c()
    geneSets = rredis::redisGet('geneSets')
    for (i in c(1:length(geneSets))) {
        pilot_genes <- union(pilot_genes, geneSets[[i]][['value']])
    }
    message(paste("Loading a pilot set of", length(pilot_genes), "genes."))
    tryCatch({
        setup(genes=pilot_genes,write.to.redis = TRUE)
        # function calls via opencpu won't be able to re-use this connection
        # just close it
        rredis::redisClose()},
        error = function(e) {
            print(traceback())
            stop(e)
        }
    )
}
