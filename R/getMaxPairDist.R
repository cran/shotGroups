## find point pair with the maximum pairwise distance
getMaxPairDist <-
function(xy) {
    UseMethod("getMaxPairDist")
}

getMaxPairDist.data.frame <-
function(xy) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getMaxPairDist")
}

getMaxPairDist.default <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(nrow(xy) < 2) {
        warning("Maximum pairwise distance needs >= 2 points")
        return(list(d=0, idx=1))
    }

    n    <- nrow(xy)                     # number of observations
    dMat <- dist(xy, method="euclidean") # distance matrix
    idx  <- which.max(as.matrix(dMat))   # maximum distance
    i    <- (idx-1) %/% n+1              # column -> point 1
    j    <- (idx-1) %%  n+1              # row    -> point 2

    return(list(d=max(dMat), idx=c(i, j)))
}
