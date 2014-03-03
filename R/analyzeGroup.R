analyzeGroup <-
function(DF, xyTopLeft=TRUE, conversion="m2cm", bandW=0.5,
         CEPtype="CorrNormal", bootCI=c("basic", "bca")) {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    #####-----------------------------------------------------------------------
    ## make sure DF has the required variable names
    varNames <- names(DF)                # what variables are present
    needsXY1 <- c("Point.X", "Point.Y")  # coordinates must have this name
    needsXY2 <- c("X", "Y")              # or this
    wantsDst <- "Distance"               # useful
    wantsAIM <- c("Aim.X", "Aim.Y")      # useful
    hasXY1   <- needsXY1 %in% varNames
    hasXY2   <- needsXY2 %in% toupper(varNames)
    hasDst   <- wantsDst %in% varNames   # useful ones we have
    hasAIM   <- wantsAIM %in% varNames   # useful ones we have

    if(!xor(all(hasXY1), all(hasXY2))) {
        stop("Coordinates must be named either X, Y or Point.X, Point.Y")
    }

    if(!all(hasDst)) {
        warning(c("The data frame is missing variable\n",
                  paste(wantsDst[!hasDst], collapse=" "), "\n",
                  "Distance is assumed to be 100"))
        DF$Distance <- 100
    }

    if(!all(hasAIM)) {
        warning(c("The data frame is missing variable(s)\n",
            paste(wantsAIM[!hasAIM], collapse=" "), "\n",
            "Point of Aim is assumed to be in (0,0)"))
        DF$Aim.X <- 0
        DF$Aim.Y <- 0
    }

    #####-----------------------------------------------------------------------
    ## prepare data: get (x,y)-coords relative to point of aim as matrix
    xy <- getXYmat(DF, xyTopLeft=xyTopLeft)

    #####-----------------------------------------------------------------------
    ## assess shape, location and spread
    dstTrgt  <- mean(DF$Distance)       # distance to target
    shape    <- groupShape(xy, plots=TRUE, bandW=bandW, outlier="mcd",
                           dstTarget=dstTrgt, conversion=conversion)
    location <- groupLocation(xy, plots=FALSE, bootCI=bootCI, level=0.95,
                              dstTarget=dstTrgt, conversion=conversion)
    spread   <- groupSpread(xy, plots=TRUE, level=0.95, CEPtype=CEPtype,
                            bootCI=bootCI, dstTarget=dstTrgt, conversion=conversion)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(c(shape, location, spread))
}
