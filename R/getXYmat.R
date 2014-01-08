getXYmat <-
function(DF, xyTopLeft=TRUE, relPOA=TRUE) {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    ## make sure DF has the required variable names
    ## coordinates need to be called either X, Y order Point.X, Point.Y
    varNames <- names(DF)                # what variables are present
    needsXY1 <- c("Point.X", "Point.Y")  # coordinates must have this name
    needsXY2 <- c("X", "Y")              # or this
    wantsAIM <- c("Aim.X", "Aim.Y")      # useful
    hasXY1   <- needsXY1 %in% varNames
    hasXY2   <- needsXY2 %in% toupper(varNames)
    hasAIM   <- wantsAIM %in% varNames   # useful ones we have

    if(!xor(all(hasXY1), all(hasXY2))) {
       stop("Coordinates must be named either X, Y or Point.X, Point.Y\n")
    }

    ## analysis should be relative to POA, but POA is missing
    if(!all(hasAIM) & relPOA) {
        warning(c("\nThe data frame is missing variable(s)\n",
                  paste(wantsAIM[!hasAIM], collapse=" "), "\n",
                  "Point of Aim is assumed to be in (0,0)\n"))
        relPOA <- FALSE
    }

    if(!relPOA) {                        # coords not relative to POA
        DF$Aim.X <- 0                    # -> set POA to (0,0)
        DF$Aim.Y <- 0  
    }
    
    ## if names are X, Y rename to Point.X, Point.Y
    if(all(hasXY2)) {
        dfNames <- names(DF)
        dfNames[dfNames %in% c("x", "X")] <- "Point.X"
        dfNames[dfNames %in% c("y", "Y")] <- "Point.Y"
        names(DF) <- dfNames
        warning("Variables X, Y were renamed to Point.X, Point.Y\n")
    }
    
    ## coords relative to point of aim
    ## y-coords exported from OnTarget: (0,0) is top-left
    X <- DF$Point.X - DF$Aim.X           # x-coords
    if(xyTopLeft) {
        Y <- -(DF$Point.Y - DF$Aim.Y)
    } else {
        Y <-   DF$Point.Y - DF$Aim.Y
    }

    return(cbind(X, Y))                  # new (x,y)-coords as matrix
}
