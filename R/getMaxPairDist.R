getMaxPairDist <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(nrow(xy) < 2) {
        warning("max pairwise distance needs >= 2 points")
        return(list(d=0, idx=1))
    }

    dMat <- dist(xy, method="euclidean")   # distance matrix
    n    <- nrow(xy)                       # number of observations
    rowN <- (n-1):1                        # number of rows in dMat

    ## find point pair with the maximum pairwise distance
    idx <- which.max(as.matrix(dMat))[1]   # maximum distance
    i   <- (idx-1) %/% n+1                 # column -> point 1
    j   <- (idx-1) %%  n+1                 # row    -> point 2

    return(list(d=max(dMat), idx=c(i, j)))
}
