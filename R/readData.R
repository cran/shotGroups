getFileNames <-
function(fPath=getwd(), fNames, fPat) {
    files <- character(0)
    ## do we have file names or a name pattern?
    if(!missing(fNames)) {               # we have file names
        files <- paste(fPath, fNames, sep="/")
    } else if(!missing(fPat)) {          # we have a name pattern
        files <- list.files(path=fPath, pattern=fPat, full.names=TRUE)
    ## are we are in interactive mode and under Windows?
    } else if(interactive() & (.Platform$OS.type == "windows")) {
        ## choose files interactively
        myFilt <- rbind(Filters, txtCsvDat=c("Data files (*.txt, *.csv, *.dat)",
                                             "*.txt;*.csv;*.dat"))
        files <- choose.files(filters=myFilt[c("txtCsvDat", "All"), ], index=1)
    }

    return(files)
}

readDataMisc <-
function(fPath=getwd(), fNames, fPat, combine=FALSE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1) { stop("no files were selected") }

    ## files should be whitespace-delimited, variable names should not have spaces
    ## read in files into a list of data frames
    DFs <- lapply(files, function(f) { read.table(f, header=TRUE) } )
    nDF <- length(DFs)                   # number of list components
    names(DFs) <- paste("file", seq(along=numeric(nDF)), sep="")  # name them

    ## make sure that the data frames all have the correct variables
    wants    <- c("Group", "Distance", "Aim.X", "Aim.Y", "Point.X", "Point.Y")
    varNames <- unlist(sapply(DFs, names))   # all variable names

    ##  build shared set of variable names
    for(i in seq(along=numeric(nDF))) {
        varNames <- intersect(varNames, names(DFs[[i]]))
    }

    has <- wants %in% varNames
    if(!all(has)) {
        warning(cat("at least one file is missing variable(s)", wants[!has],
                    "\nthat may be required later by analysis functions\n"))
    }

    if(combine) { return(combineData(DFs)) } else { return(DFs) }
}

readDataOT1 <-
function(fPath=getwd(), fNames, fPat, combine=FALSE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1) { stop("no files were selected") }

    ## files should be original OnTarget-output (version 1.**) with variables
    ## Project Title, Group, Ammunition, Distance, Aim X, Aim Y,
    ## Center X, Center Y, Point X, Point Y (10 fields + trailing tab = 11)
    nFields <- unlist(lapply(files, function(x) count.fields(x, sep="\t")))
    if(!all(nFields == 11)) {
        stop(
"it appears at least one file does not contain exactly
the required set of 10 variables - see help(readDataOT1)
maybe you should use readDataMisc instead")
    }

    ## read in files into a list of data frames
    DFs <- lapply(files, function(f) {
              read.delim(f, colClasses=c("character", "factor", "character",
                         "numeric", "numeric", "numeric", "numeric", "numeric",
                         "numeric", "numeric", "NULL"), strip.white=TRUE) } )
    nDF <- length(DFs)                   # number of list components
    names(DFs) <- paste("file", seq(along=numeric(nDF)), sep="")  # name them

    ## make sure that the data frames all have the correct variables
    wants    <- c("Group", "Distance", "Aim.X", "Aim.Y", "Point.X", "Point.Y")
    varNames <- unlist(sapply(DFs, names))   # all variable names

    ##  build shared set of variable names
    for(i in seq(along=numeric(nDF))) {
        varNames <- intersect(varNames, names(DFs[[i]]))
    }

    has <- wants %in% varNames
    if(!all(has)) {
        warning(cat("at least one file is missing variable(s)", wants[!has],
                    "\nthat may be required later by analysis functions\n"))
    }

    if(combine) { return(combineData(DFs)) } else { return(DFs) }
}

readDataOT2 <-
function(fPath=getwd(), fNames, fPat, combine=FALSE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1) { stop("no files were selected") }

    ## files should be original OnTarget-output (version 2.**) with variables
    ## Project Title, Group, Ammunition, Distance, Aim X, Aim Y,
    ## Center X, Center Y, Point X, Point Y (10 fields)
    nFields <- unlist(lapply(files, function(x) count.fields(x, sep=",")))
    if(!all(nFields == 10)) {
        stop(
"it appears at least one file does not contain exactly
the required set of 10 variables - see help(readDataOT2)
maybe you should use readDataMisc instead")
    }

    ## read in files into a list of data frames
    DFs <- lapply(files, function(f) {
              read.csv(  f, colClasses=c("character", "factor", "character",
                         "numeric", "numeric", "numeric", "numeric", "numeric",
                         "numeric", "numeric"),         strip.white=TRUE) } )
    nDF <- length(DFs)                   # number of list components
    names(DFs) <- paste("file", seq(along=numeric(nDF)), sep="")  # name them

    ## make sure that the data frames all have the correct variables
    wants    <- c("Group", "Distance", "Aim.X", "Aim.Y", "Point.X", "Point.Y")
    varNames <- unlist(sapply(DFs, names))   # all variable names

    ##  build shared set of variable names
    for(i in seq(along=numeric(nDF))) {
        varNames <- intersect(varNames, names(DFs[[i]]))
    }

    has <- wants %in% varNames
    if(!all(has)) {
        warning(cat("at least one file is missing variable(s)", wants[!has],
                    "\nthat may be required later by analysis functions\n"))
    }

    if(combine) { return(combineData(DFs)) } else { return(DFs) }
}
