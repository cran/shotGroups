getMinBBox <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }
    if(nrow(xy) < 2)    { stop("xy must have at least two rows") }

    ## rotating pliers algorithm using the convex hull
    H    <- chull(xy)                    # hull indices, vertices ordered clockwise
    n    <- length(H)                    # number of hull vertices
    hull <- xy[H, ]                      # hull vertices

    ## unit basis vectors for all subspaces spanned by the hull edges
    hDir  <- diff(rbind(hull, hull[1,])) # account for circular hull vertices
    hLens <- sqrt(rowSums(hDir^2))       # length of basis vectors
    huDir <- diag(1/hLens) %*% hDir      # scaled to unit length

    ## unit basis vectors for the orthogonal subspaces
    ## rotation by 90 deg -> y' = x, x' = -y
    ouDir <- cbind(-huDir[ , 2], huDir[ , 1])

    ## project hull vertices on the subspaces spanned by the hull edges, and on
    ## the subspaces spanned by their orthogonal complements - in subspace coords
    projHull <- huDir %*% t(hull)    # row = hull edge subspace, col = hull vert
    projOrth <- ouDir %*% t(hull)    # row = orth subspace,      col = hull vert

    ## extreme projections that mark rectangle corners
    extrIdxH <- t(apply(projHull, 1, function(x) { c(which.min(x), which.max(x)) } ))
    extrIdxO <- t(apply(projOrth, 1, function(x) { c(which.min(x), which.max(x)) } ))

    ## areas of bounding rectangles (in subspace coords)
    widths  <- abs(diff(rbind(projHull[cbind(1:n, extrIdxH[ , 1])],
                              projHull[cbind(1:n, extrIdxH[ , 2])])))
    heights <- abs(diff(rbind(projOrth[cbind(1:n, extrIdxO[ , 1])],
                              projOrth[cbind(1:n, extrIdxO[ , 2])])))

    ## hull edge leading to the minimum-area bounding rectangle
    eMin <- which.min(widths*heights)
    w    <- widths[eMin]                 # width of minimum bounding rect
    h    <- heights[eMin]                # height of minimum bounding rect

    ## extreme projections in subspace coordinates
    hProj <- rbind(   projHull[eMin, extrIdxH[eMin, ]], 0)
    oProj <- rbind(0, projOrth[eMin, extrIdxO[eMin, ]])

    ## move projections to rectangle corners
    hPts <- sweep(hProj, 1, oProj[ , 1], "+")
    oPts <- sweep(hProj, 1, oProj[ , 2], "+")

    ## corners in standard coordintes
    basis <- cbind(huDir[eMin, ], ouDir[eMin, ])  # basis formed by hull edge and orth
    hCorn <- basis %*% hPts
    oCorn <- basis %*% oPts

    ## angle of longer edge pointing up
    e   <- if(w > h)  { hCorn[, 2] - hCorn[, 1] } else { oCorn[, 2] - oCorn[, 1] }
    eUp <- e * sign(e[2])                  # rotate upwards 180 deg if necessary
    deg <- atan2(eUp[2], eUp[1])*180 / pi  # angle in degrees

    ## return (4x2)-matrix: reverse point order -> make it usable in polygon()
    pts <- t(cbind(hCorn, oCorn[ , c(2, 1)]))

    return(list(pts=pts, width=w, height=h, angle=deg))
}
