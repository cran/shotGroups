getMinBBox <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }
    if(nrow(xy) < 2)    { stop("xy must have at least two rows") }

    ## algorithm using the convex hull
    H    <- chull(xy)                    # convex hull
    n    <- length(H)                    # number of hull vertices
    Hidx <- seq(along=numeric(n))        # index for vertices
    post <- (Hidx %% n) + 1              # next vertex in S
    area <- numeric(n)                   # for rectangle areas

    ## some objects we later need to collect information
    VV   <- array(numeric(n*n*2), c(n, n, 2))
    UU   <- array(numeric(n*n*2), c(n, n, 2))
    V    <- matrix(numeric(n*n), ncol=n)
    U    <- matrix(numeric(n*n), ncol=n)
    vIdx <- matrix(numeric(n*2), ncol=2)
    uIdx <- matrix(numeric(n*2), ncol=2)
    lenA <- numeric(n)     # for lengths of basis vectors of subspaces
    lenB <- numeric(n)     # for lengths of basis vectors of subspaces

    ## iterate over all hull edges
    for(i in Hidx) {
        x <- xy[H[i], ]
        y <- xy[H[post[i]], ]

        ## subspace spanned by hull edge and its orthogonal complement
        a       <- y - x                 # hull edge = basis vector
        lenA[i] <- sqrt(sum(a^2))        # its length
        v       <- a %*% solve(t(a) %*% a) %*% t(a) %*% x  # projection onto edge
        b       <- x-v                   # spans subspace orthogonal to hull edge
        lenB[i] <- sqrt(sum(b^2))        # length basis vector

        ## project all vertices on these subspaces
        for(j in Hidx) {
            z          <- xy[H[j], ]     # current hull vertex
            V[i, j]    <- solve(t(a) %*% a) %*% t(a) %*% z  # subspace coords
            U[i, j]    <- solve(t(b) %*% b) %*% t(b) %*% z  # orth subspace coords
            VV[i, j, ] <- a * V[i, j]    # coords in original system
            UU[i, j, ] <- b * U[i, j]    # coords in original system
        }

        ## store vertices that give extreme projections, defining the rectangle
        vIdx[i, ] <- c(which.min(V[i, ]), which.max(V[i, ]))
        uIdx[i, ] <- c(which.min(U[i, ]), which.max(U[i, ]))

        ## store area of bounding rectangle (in subspace coords)
        area[i] <- abs(diff(V[i, vIdx[i, ]])) *
                   abs(diff(U[i, uIdx[i, ]]))
    }

    ## hull vertex leading to the minimum-area rectangle
    iMin <- which.min(area)

    ## width and height from extreme subspace projections
    ## in original coordinates
    w <- abs(diff(V[iMin, vIdx[iMin, ]])) * lenA[iMin]
    h <- abs(diff(U[iMin, uIdx[iMin, ]])) * lenB[iMin]

    ## extreme subspace projections in original coordinates
    vMat <- VV[iMin, vIdx[iMin, ], ]
    uMat <- UU[iMin, uIdx[iMin, ], ]

    ## angle of longer edge pointing up
    e   <- if(w > h)    { vMat[2, ] - vMat[1, ] } else { uMat[2, ] - uMat[1, ] }
    e2  <- if(e[2] < 0) { -e } else { e }
    deg <- atan2(e2[2], e2[1])*180 / pi  # angle in degree

    ## move projections to envelope point set
    p1 <- sweep(vMat, 2, uMat[1, ], "+")
    p2 <- sweep(vMat, 2, uMat[2, ], "+")

    ## returned matrix: reverse point order -> usable in polygon()
    pts <- rbind(p1, p2[c(2, 1), ])
    # xy[H[vIdx[iMin, ]], ]
    # xy[H[uIdx[iMin, ]], ]

    return(list(pts=pts, width=w, height=h, angle=deg))
}
