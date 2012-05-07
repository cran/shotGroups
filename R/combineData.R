combineData <-
function(DFs) {
    if(!is.list(DFs)) { stop("DFs must be a list") }
    if(!all(sapply(DFs, is.data.frame))) { stop("DFs must be a list of data frames") }

    ## restrict data frames to those variables that are shared by them all
    nDF      <- length(DFs)                  # number of data frames in list
    varNames <- unlist(sapply(DFs, names))   # all variable names
    for(i in seq(along=numeric(nDF))) {      # get shared set
        varNames <- intersect(varNames, names(DFs[[i]]))
    }

    ## check if data frames contain required variables
    needs <- c("Group")                      # absolutely required
    wants <- c("Distance", "Aim.X", "Aim.Y", "Point.X", "Point.Y")  # useful
    has1  <- needs %in% varNames             # required ones we have
    has2  <- wants %in% varNames             # useful ones we have
    if(!all(has1)) {
        stop(cat("at least one data frame is missing variable(s)\n", needs[!has1], "\n"))
    }

    if(!all(has2)) {
        warning(cat("at least one data frame is missing variable(s)\n", wants[!has2],
                    "\nthat may be required later by analysis functions\n"))
    }

    DFrestr <- lapply(DFs, function(x) subset(x, select=varNames))
    nObs    <- sapply(DFrestr, nrow)         # number of observations in each data frame
    DFall   <- do.call("rbind", DFrestr)     # combine data frames
    rownames(DFall) <- NULL                  # remove row names

    ## add new factor Origin for coding original file
    DFall$Origin <- factor(rep(seq(along=numeric(nDF)), nObs))

    ## add new factor Series for coding Groups as a consecutive number over files
    ## first a factor with alphabetically ordered levels
    DFall$orgSer <- droplevels(interaction(DFall$Origin, DFall$Group))

    ## convert orgSer to a factor that with consecutively numbered levels
    runs         <- rle(as.character(DFall$orgSer))
    runs$values  <- 1:length(runs$values)
    DFall$Series <- factor(inverse.rle(runs))

    return(DFall)
}
