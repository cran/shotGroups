getMaxPairDist <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }

    dMat <- dist(xy, method="euclidean")   # distance matrix
    n    <- nrow(xy)                       # number of observations
    rowN <- (n-1):1                        # number of rows in dMat

    ## find point pair with the maximum pairwise distance
    ind <- which.max(dMat)                 # linear index
    i   <- sum(cumsum(rowN) < ind) + 1     # first point
    j   <- n - (sum(rowN[1:i]) - ind)      # second point

    return(list(d=max(dMat), idx=c(i, j)))
}
