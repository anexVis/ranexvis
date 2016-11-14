.onLoad <- function(libname, pkgname) {
    # open redis connection
    rredis::redisConnect('localhost')

    ### for development only
    # loading full set of expression data is not possible with the current implementation and setup:
    # redis server has a limit of 2GB-value size, and redis cluster has to be installed on multiple physical node
    pilot_genes = c('')
    # setup(write.to.redis = TRUE)
}
