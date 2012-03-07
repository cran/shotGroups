getDistToCtr <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }

    xyCtr <- scale(xy, scale=FALSE, center=TRUE)  # centered data
    return(sqrt(rowSums(xyCtr^2)))                # distances to center
}
