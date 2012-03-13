getConfEll <-
function(xy, level=0.5, dstTarget=25, conversion="m2cm") {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }
    if(!is.numeric(level))          { stop("level must be numeric") }
    if((level <= 0) | (level >= 1)) { stop("level must be in (0,1)") }

    ## check for minimum number of points
    if(nrow(xy) < 2) { stop("we need >= 2 points for confidence ellipse") }

    haveRob <- TRUE                      # can we do robust estimation?
    if(nrow(xy) < 4) {
        warning("we need >= 4 points for robust estimations")
        haveRob <- FALSE
    }

    ## group center and covariance matrix
    covXY  <- cov(xy)                    # covariance matrix (x,y)-coords
    ellRad <- sqrt(eigen(covXY)$values)  # radii error ellipse

    ## magnification factor error ellipse -> confidence ellipse
    N    <- nrow(xy)                     # number of observations
    dfn  <- ncol(xy)                     # numerator df
    dfd  <- N-1                          # denominator df
    mag  <- sqrt(dfn*qf(level, dfn, dfd))  # magnification factor = t-value
    size <- rbind(unit=mag*ellRad,       # radii confidence ellipse
                  MOA=getMOA(mag*ellRad, dstTarget, conversion))
    colnames(size) <- c("semi-major", "semi-minor")

    ## aspect ratio and flattening of characteristic ellipse
    aspRatio <- ellRad[1] / ellRad[2]
    flat     <- 1 - (ellRad[2] / ellRad[1])
    shape    <- c(aspectRatio=aspRatio, flattening=flat)

    if(haveRob) {                        # same for robust estimation
        # library(robustbase)                  # for covMcd()
        rob       <- covMcd(xy)
        covXYrob  <- rob$cov             # robust estimate: group covariance matrix
        ellRadRob <- sqrt(eigen(covXYrob)$values)  # radii error ellipse
        sizeRob   <- rbind(unit=mag*ellRadRob,     # radii robust confidence ellipse
                           MOA=getMOA(mag*ellRadRob, dstTarget, conversion))
        colnames(sizeRob) <- c("semi-major", "semi-minor")

        aspRatioRob <- ellRadRob[1] / ellRadRob[2]
        flatRob     <- 1 - (ellRadRob[2] / ellRadRob[1])
        shapeRob    <- c(aspectRatio=aspRatioRob, flattening=flatRob)
    } else { sizeRob=NULL; shapeRob=NULL }

    return(list(size=size,   sizeRob=sizeRob,
                shape=shape, shapeRob=shapeRob, magFac=mag))
}
