## circle defined by three points
getCircleFrom3 <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if((nrow(xy) != 3) | (ncol(xy) != 2)) { stop("xy must be (3x2)-matrix") }

    aa <- xy[1,  ]
    bb <- xy[2,  ]
    cc <- xy[3,  ]
    y  <- xy[ , 2]

    xDeltaA <- bb[1] - aa[1]
    yDeltaA <- bb[2] - aa[2]
    xDeltaB <- cc[1] - bb[1]
    yDeltaB <- cc[2] - bb[2]
    xDeltaC <- cc[1] - aa[1]
    yDeltaC <- cc[2] - aa[2]

    ## check if the points are collinear: qr(xy)$rank == 1, or:
    ## determinant of difference matrix = 0, no need to use det()
    dMat <- rbind(c(xDeltaA, yDeltaA), c(xDeltaB, yDeltaB))
    if(isTRUE(all.equal(dMat[1,1]*dMat[2,2] - dMat[1,2]*dMat[2,1], 0, check.attributes=FALSE))) {
        ## define the circle as the one that's centered between them
        rangeX <- range(c(aa[1], bb[1], cc[1]))
        rangeY <- range(c(aa[2], bb[2], cc[2]))
        ctr    <- c(rangeX[1] + 0.5*diff(rangeX), rangeY[1] + 0.5*diff(rangeY))
        rad    <- sqrt((0.5*diff(rangeX))^2 + (0.5*diff(rangeY))^2)
    } else {
        rad <- prod(dist(xy)) / (2 * abs(det(cbind(xy, 1))))  # circle radius
        v1  <- rowSums(xy^2)                    # first vector in the numerator
        v2x <- c( xDeltaB, -xDeltaC,  xDeltaA)  # 2nd vector numerator for Mx
        v2y <- c(-yDeltaB,  yDeltaC, -yDeltaA)  # 2nd vector numerator for My
        ctr <- c(t(v1) %*% v2y, t(v1) %*% v2x) / (2 * (t(y) %*% v2x))  # center
    }

    return(list(ctr=ctr, rad=rad))
}

## alternative formulas
## detABC <- (cx-ax)*(cy+ay) + (bx-cx)*(by+cy) + (ax-bx)*(ay+by)
## Mx <- 0.5 * (((ax^2 + ay^2)*(by-cy) + (bx^2+by^2)*(cy-ay) + (cx^2+cy^2)*(ay-by)) /
##              (ay*(cx-bx) + by*(ax-cx) + cy*(bx-ax)))
## My <- 0.5 * (((ax^2 + ay^2)*(cx-bx) + (bx^2+by^2)*(ax-cx) + (cx^2+cy^2)*(bx-ax)) /
##              (ay*(cx-bx) + by*(ax-cx) + cy*(bx-ax)))
## rad <- dist(rbind(xy, ctr))[3]

## vertex that produces the circles with the maximum radius
## used in getMinCircle()
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

## angle at B in triangle ABC
getAngleTri <- function(xy, deg=TRUE) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if((nrow(xy) != 3) | (ncol(xy) != 2)) { stop("xy must be (3x2)-matrix") }

    d   <- dist(xy)
    dAB <- d[1]
    dAC <- d[2]
    dBC <- d[3]

    ## dAB*dAC should not be 0
    if(isTRUE(all.equal(dAB*dAC, 0, check.attributes=FALSE))) {
        stop("some edges have zero length")
    }

    Wabc <- (dAB^2 + dBC^2 - dAC^2)
    arc  <- acos(Wabc / (2*dAB*dBC))      # angle in radians
    ang  <- ifelse(deg, arc*180/pi, arc)  # return angle in degree or radians

    return(ang)
}

## checks if the angle at B in triangle ABC is bigger than 90 degrees
isBiggerThan90 <- function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if((nrow(xy) != 3) | (ncol(xy) != 2)) { stop("xy must be (3x2)-matrix") }

    d   <- dist(xy)
    dAB <- d[1]
    dAC <- d[2]
    dBC <- d[3]
    return((dAB^2 + dBC^2 - dAC^2) < 0)
}

## minimal enclosing circle
getMinCircle <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(nrow(xy) < 2)    { stop("there must be at least two points") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }

    ## Skyum algorithm based on the convex hull
    H <- chull(xy)                       # hull indices (vertices ordered clockwise)
    S <- H                               # copy that will be changed
    while(length(S) >= 2) {
        n    <- length(S)                # number of remaining hull vertices
        Sidx <- seq(along=numeric(n))    # index for vertices
        post <- (Sidx %% n) + 1          # next vertex in S
        prev <- Sidx[order(post)]        # previous vertex in S
        mIdx <- getMaxRad(xy, S)         # idx for maximum radius

        ## triangle where mIdx is vertex B in ABC
        Smax <- rbind(xy[S[prev[mIdx]], ], xy[S[mIdx], ], xy[S[post[mIdx]], ])

        ## if there's only two hull vertices, we're done
        if(n <= 2) { break }

        ## check if angle(ABC) is > 90
        ## if so, eliminate B - if not, we're done
        if(isBiggerThan90(Smax)) { S <- S[-mIdx] } else { break }
    }

    return(getCircleFrom3(Smax))
}
