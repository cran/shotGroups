analyzeGroup <-
function(DF, xyTopLeft=TRUE, conversion="m2cm", bandW=0.5) {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    #####-----------------------------------------------------------------------
    ## make sure DF has the required variable names
    varNames <- names(DF)                # what variables are present
    needs    <- c("Distance", "Aim.X", "Aim.Y", "Point.X", "Point.Y")   # what we need
    has      <- needs %in% varNames      # what we have of the required ones
    if(!all(has)) {
        stop(cat("the data frame is missing variable(s)\n", needs[!has], "\n"))
    }

    #####-----------------------------------------------------------------------
    ## prepare data: coords relative to point of aim
    ## y-coords exported from OnTarget: (0,0) is top-left
    X <- DF$Point.X - DF$Aim.X           # x-coords
    if(xyTopLeft) {
        Y <- -(DF$Point.Y - DF$Aim.Y)
    } else {
        Y <-   DF$Point.Y - DF$Aim.Y
    }
    xy <- cbind(X, Y)                    # new (x,y)-coords as matrix

    #####-----------------------------------------------------------------------
    ## assess shape, location and spread
    dstTrgt  <- mean(DF$Distance)       # distance to target
    shape    <- groupShape(   xy, plots=TRUE, bandW=bandW)
    location <- groupLocation(xy, plots=FALSE, dstTarget=dstTrgt, conversion=conversion)
    spread   <- groupSpread(  xy, plots=TRUE, level=0.5, dstTarget=dstTrgt,
                            conversion=conversion)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(c(shape, location, spread))
}
