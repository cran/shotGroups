getXYmat <-
function(DF, xyTopLeft=TRUE, relPOA=TRUE) {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    ## make sure DF has the required variable names
    ## coordinates need to be called either X, Y order Point.X, Point.Y
    dfNames  <- names(DF)                # what variables are present
    needsXY1 <- c("Point.X", "Point.Y")  # coordinates must have this name
    needsXY2 <- c("X", "Y")              # or this
    wantsAIM <- c("Aim.X", "Aim.Y")      # useful
    hasXY1   <- needsXY1 %in% dfNames
    hasXY2   <- needsXY2 %in% toupper(dfNames)
    hasAIM   <- wantsAIM %in% dfNames   # useful ones we have

    if(!xor(all(hasXY1), all(hasXY2))) {
       stop("Coordinates must be named either X, Y or Point.X, Point.Y")
    }

    if(("Z" %in% toupper(dfNames)) && ("Point.Z" %in% dfNames)) {
        stop("Coordinates must be named either Z or Point.Z")
    }

    ## analysis should be relative to POA, but POA is missing
    if(!all(hasAIM) && relPOA) {
        warning(c("The data frame is missing variable(s)\n",
                  paste(wantsAIM[!hasAIM], collapse=" "), "\n",
                  "Point of Aim is assumed to be in (0,0)"))
        relPOA <- FALSE
    }

    if(!relPOA) {                        # coords not relative to POA
        DF$Aim.X <- 0                    # -> set POA to (0,0)
        DF$Aim.Y <- 0

        if(("Z" %in% toupper(dfNames)) || ("Point.Z" %in% dfNames)) {
            DF$Aim.Z <- 0
        }
    }
    
    ## if names are X, Y rename to Point.X, Point.Y
    if(all(hasXY2)) {
        dfNames <- names(DF)
        dfNames[dfNames %in% c("x", "X")] <- "Point.X"
        dfNames[dfNames %in% c("y", "Y")] <- "Point.Y"
        dfNames[dfNames %in% c("z", "Z")] <- "Point.Z"
        names(DF) <- dfNames
        warning("Variables X, Y were renamed to Point.X, Point.Y")
    }
    
    ## coords relative to point of aim
    ## y-coords exported from OnTarget: (0,0) is top-left
    X <- DF$Point.X - DF$Aim.X           # x-coords
    if(xyTopLeft) {
        Y <- -(DF$Point.Y - DF$Aim.Y)
    } else {
        Y <-   DF$Point.Y - DF$Aim.Y
    }

    if(("Z" %in% toupper(dfNames)) || ("Point.Z" %in% dfNames)) {
        Z <- DF$Point.Z - DF$Aim.Z       # x-coords
    } else {
        Z <- NULL
    }

    return(cbind(X, Y, Z))               # new (x,y)-coords as matrix
}
