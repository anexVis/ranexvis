.onLoad <- function(libname, pkgname) {
    # open redis connection
    rredis::redisConnect('localhost')
    setup(write.to.redis = TRUE)
}
