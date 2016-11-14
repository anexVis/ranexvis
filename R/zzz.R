.onLoad <- function(libname, pkgname) {
    # open redis connection
    # nodelay=TRUE cause Error in setting TCP_NODELAY
    rredis::redisConnect('localhost',nodelay = FALSE)

    ### for development only
    # loading full set of expression data is not possible with the current implementation and setup:
    # redis server has a limit of 2GB-value size, and redis cluster has to be installed on multiple physical node
    pilot_genes = c('XYLT1', 'XYLT2','B4GALT7', 'B3GALT6',
                    'B3GAT3', 'EXTL2', 'EXTL3', 'EXTL1',
                    'EXT1', 'EXT2', 'NDST1', 'NDST2',
                    'NDST3', 'NDST4', 'GLCE', 'HS2ST1',
                    'HS6ST1', 'HS6ST2', 'HS6ST3', 'HS3ST1',
                    'HS3ST2', 'HS3ST3A1', 'HS3ST3B1', 'HS3ST4',
                    'HS3ST5', 'HS3ST6')
    tryCatch(
        setup(genes=hgnc2ensembl(pilot_genes),write.to.redis = TRUE),
        error = function(e) {
            print(traceback())
            stop(e)
        }
    )
}
