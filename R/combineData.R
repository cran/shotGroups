combineData <-
function(DFs) {
    if(!is.list(DFs)) { stop("DFs must be a list") }
    if(!all(sapply(DFs, is.data.frame))) { stop("DFs must be a list of data frames") }

    ## shared set of variable names
    varNameL <- lapply(DFs, names)       # list of variable names
    varNames <- Reduce(intersect, varNameL)

    ## check if data frames contain required variables
    wantsGrp <- "Group"                  # useful
    wantsDst <- "Distance"               # useful
    hasGrp   <- wantsGrp %in% varNames   # useful ones we have
    hasDst   <- wantsDst %in% varNames   # useful ones we have

    if(!all(hasGrp)) {
        warning(c("At least one data frame is missing variable\n",
                  paste(wantsGrp[!hasGrp], collapse=" "),
                  "\nGroup is set to 1\n"))
    
        setGroup <- function(x) {
            if(!("Group" %in% names(x))) { x$Group <- 1 }
            x
        }
        DFs <- lapply(DFs, setGroup)
    }

    if(!all(hasDst)) {
        warning(c("At least one file is missing variable\n",
                  paste(wantsDst[!hasDst], collapse=" "),
                  "\nthat may be required later by analysis functions\n"))
    }

    ## make sure each data frame has either X, Y or Point.X, Point.Y
    replaceXY <- function(x) {
        dfNames  <- names(x)
        needsXY1 <- c("Point.X", "Point.Y")  # coordinates must have this name
        needsXY2 <- c("X", "Y")              # or this
        hasXY1   <- needsXY1 %in% dfNames
        hasXY2   <- needsXY2 %in% toupper(dfNames)
        
        if(!xor(all(hasXY1), all(hasXY2))) { # not (either X, Y or Point.X, Point.Y)
            stop("Coordinates must be named either X, Y or Point.X, Point.Y\n")
        }
        
        ## if X, Y -> rename to Point.X, Point.Y
        if(all(hasXY2)) {
            dfNames[dfNames %in% c("x", "X")] <- "Point.X"
            dfNames[dfNames %in% c("y", "Y")] <- "Point.Y"
            names(x) <- dfNames
            warning("Variables X, Y were renamed to Point.X, Point.Y\n")
        }
        x
    }
    
    DFs <- lapply(DFs, replaceXY)
    
    ## restrict data frames to shared variables variables
    varsNow <- Reduce(intersect, lapply(DFs, names))  # shared set of variables
    DFrestr <- lapply(DFs, function(x) x[, varsNow])  # only select these
    nObs    <- sapply(DFrestr, nrow)         # number of observations in each data frame
    DFall   <- do.call("rbind", DFrestr)     # combine data frames
    rownames(DFall) <- NULL                  # remove row names

    ## add new factor Origin for coding original file
    DFall$Origin <- factor(rep(seq(along=DFs), nObs))

    ## add new factor Series for coding Groups as a consecutive number over files
    ## first a factor with alphabetically ordered levels
    DFall$orgSer <- droplevels(interaction(DFall$Origin, DFall$Group))

    ## convert orgSer to a factor with consecutively numbered levels
    runs         <- rle(as.character(DFall$orgSer))
    runs$values  <- 1:length(runs$values)
    DFall$Series <- factor(inverse.rle(runs))

    return(DFall)
}
