getXYmat <-
function(DF, xyTopLeft=TRUE) {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    ## make sure DF has the required variable names
    varNames <- names(DF)                # what variables are present
    needs    <- c("Aim.X", "Aim.Y", "Point.X", "Point.Y")   # what we need
    has      <- needs %in% varNames      # what we have of the required ones
    if(!all(has)) {
        stop(cat("the data frame is missing variable(s)\n", needs[!has], "\n"))
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
