getBoundingBox <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }

    x   <- range(xy[ , 1])
    y   <- range(xy[ , 2])
    pts <- c(xleft=x[1], ybottom=y[1], xright=x[2], ytop=y[2])

    return(list(pts=pts, width=abs(diff(x)), height=abs(diff(y))))
}
