getMaxRad <- function(xy, S) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }
    if(!is.numeric(S))  { stop("S must be numeric") }
    if((length(S) < 2) | (nrow(xy) < 2)) { stop("there must be at least two points") }
    if(length(S) > nrow(xy)) { stop("there can only be as many indices in S as points in xy") }

    n    <- length(S)                    # number of points
    Sidx <- seq(along=numeric(n))        # index for points
    rads <- numeric(n)                   # radii for all circles
    post <- (Sidx %% n) + 1              # next point in S
    prev <- Sidx[order(post)]            # previous point in S
    for(i in Sidx) {
        pts     <- rbind(xy[S[prev[i]], ], xy[S[i], ], xy[S[post[i]], ])
        rads[i] <- getCircleFrom3(pts)$rad  # circle radius
    }

    return(which.max(rads))
}
