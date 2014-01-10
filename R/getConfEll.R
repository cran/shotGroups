getConfEll <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", doRob=TRUE) {
    UseMethod("getConfEll")
}

getConfEll.data.frame <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", doRob=TRUE) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getConfEll")
}

getConfEll.default <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", doRob=TRUE) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }
    if(nrow(xy) < 2)    { stop("we need >= 2 points for confidence ellipse") }
    if(!is.numeric(level))          { stop("level must be numeric") }
    if((level <= 0) | (level >= 1)) { stop("level must be in (0,1)") }

    ## group center and covariance matrix
    ctr    <- colMeans(xy)               # group center
    covXY  <- cov(xy)                    # covariance matrix (x,y)-coords
    ellRad <- sqrt(eigen(covXY)$values)  # radii error ellipse
    trXY   <- sum(diag(covXY))           # trace of covariance matrix
    detXY  <- covXY[1,1]*covXY[2,2] - covXY[1,2]*covXY[2,1]  # determinant

    ## magnification factor error ellipse -> confidence ellipse
    N   <- nrow(xy)                      # number of observations
    dfn <- ncol(xy)                      # numerator df
    dfd <- N-1                           # denominator df
    mag <- sqrt(dfn*qf(level, dfn, dfd)) # magnification factor = t-value
    size <- rbind(unit=mag*ellRad,       # radii confidence ellipse
                  MOA=getMOA(mag*ellRad, dstTarget, conversion))
    colnames(size) <- c("semi-major", "semi-minor")

    ## aspect ratio and flattening of characteristic ellipse
    aspRatio <- ellRad[1] / ellRad[2]
    flat     <- 1 - (ellRad[2] / ellRad[1])
    shape    <- c(aspectRatio=aspRatio, flattening=flat, trace=trXY, det=detXY)

    if(nrow(xy) < 4) {
        haveRob <- FALSE
        if(doRob) {
            warning("We need >= 4 points for robust estimations")
        }
    } else {
        haveRob <- TRUE
    }                                    # if(nrow(xy) < 4)
    
    if(doRob & haveRob) {                # same for robust estimation
        rob       <- robustbase::covMcd(xy)
        ctrRob    <- rob$center
        covXYrob  <- rob$cov
        ellRadRob <- sqrt(eigen(covXYrob)$values)  # radii error ellipse
        trXYrob   <- sum(diag(covXYrob))           # trace of covariance matrix
        detXYrob  <- covXYrob[1,1]*covXYrob[2,2] - covXYrob[1,2]*covXYrob[2,1]  # determinant
        sizeRob   <- rbind(unit=mag*ellRadRob,     # radii robust confidence ellipse
                           MOA=getMOA(mag*ellRadRob, dstTarget, conversion))
        colnames(sizeRob) <- c("semi-major", "semi-minor")

        aspRatioRob <- ellRadRob[1] / ellRadRob[2]
        flatRob     <- 1 - (ellRadRob[2] / ellRadRob[1])
        shapeRob    <- c(aspectRatio=aspRatioRob, flattening=flatRob, trace=trXYrob, det=detXYrob)
    }                                    # if(doRob & haveRob)

    ## set robust estimates to NULL if not available
    if(!(doRob & haveRob)) {
        sizeRob  <- NULL
        shapeRob <- NULL
        ctrRob   <- NULL
        covXYrob <- NULL
    }                                    # if(!(doRob & haveRob))
    
    return(list(ctr=ctr, ctrRob=ctrRob, cov=covXY, covRob=covXYrob,
                size=size, sizeRob=sizeRob,
                shape=shape, shapeRob=shapeRob, magFac=mag))
}
