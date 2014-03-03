## read files that are comma or whitespace-delimited
## variable names must not have spaces
readDataMisc <-
function(fPath=".", fNames, fPat, combine=TRUE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1) { stop("No files were selected") }

    ## determine whether file is csv or whitespace delimited
    isCSV <- function(ext, nFields) {
        csv <- if(tolower(ext) == "csv") {        # explicit file extension
            TRUE
        } else if((length(unique(nFields)) == 1) && (nFields[1] >= 2)) {
            ## same number of comma-separated fields in all rows
            ## and at least two comma-separated fields
            TRUE
        } else {
            FALSE
        }
    }

    ## read a single file, possibly csv, possibly whitespace delimited
    readMe <- function(f) {
        pieces  <- strsplit(f, "\\.")[[1]]
        ext     <- tolower(pieces[length(pieces)]) # file extensions
        nFields <- count.fields(f, sep=",")
        readFun <- if(isCSV(ext, nFields)) { read.csv } else { read.table }
        readFun(f, header=TRUE)
    }

    ## read multiple files, some possibly csv, some possibly whitespace delimited
    DFs <- lapply(files, readMe)
    names(DFs) <- paste("file", seq(along=DFs), sep="")  # name them

    ## build shared set of variable names
    varNameL <- lapply(DFs, names)       # list of variable names
    varNames <- Reduce(intersect, varNameL)

    ## make sure that the data frames all have the correct variables
    wants <- c("Distance", "Group", "Aim.X", "Aim.Y")  # useful
    has   <- wants %in% varNames

    if(!all(has)) {
        warning(c("At least one file is missing variable(s)\n",
            paste(wants[!has], collapse=" "),
            "\nthat may be required later by analysis functions"))
    }

    ## make sure each data frame has either X, Y or Point.X, Point.Y
    replaceXY <- function(x) {
        dfNames  <- names(x)
        needsXY1 <- c("Point.X", "Point.Y")  # coordinates must have this name
        needsXY2 <- c("X", "Y")              # or this
        hasXY1   <- needsXY1 %in% dfNames
        hasXY2   <- needsXY2 %in% toupper(dfNames)

        if(!xor(all(hasXY1), all(hasXY2))) { # not (either X, Y or Point.X, Point.Y)
            stop("Coordinates must be named either X, Y or Point.X, Point.Y")
        }

        if(("Z" %in% toupper(dfNames)) && ("Point.Z" %in% dfNames)) {
            stop("Coordinates must be named either Z or Point.Z")
        }

        ## if X, Y, Z -> rename to Point.X, Point.Y, Point.Z
        if(all(hasXY2)) {
            dfNames[dfNames %in% c("x", "X")] <- "Point.X"
            dfNames[dfNames %in% c("y", "Y")] <- "Point.Y"
            dfNames[dfNames %in% c("z", "Z")] <- "Point.Z"
            names(x) <- dfNames
            warning("Variables X, Y were renamed to Point.X, Point.Y")
        }
        x
    }

    DFs <- lapply(DFs, replaceXY)

    if(combine) { return(combineData(DFs)) } else { return(DFs) }
}

## read files from OnTarget-output (version 1.*)
readDataOT1 <-
function(fPath=".", fNames, fPat, combine=TRUE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1) { stop("No files were selected") }

    ## assumed variables: Project Title, Group, Ammunition, Distance,
    ## Aim X, Aim Y, Center X, Center Y, Point X, Point Y
    ## 10 fields + trailing tab = 11
    nFields <- unlist(lapply(files, function(x) count.fields(x, sep="\t")))
    if(!all(nFields == 11)) {
        stop(c("It appears at least one file does not contain exactly\n",
               "the required set of 10 variables - see help(readDataOT1)\n",
               "maybe you should use readDataMisc() instead"))
    }

    ## read in files into a list of data frames
    DFs <- lapply(files, function(f) {
              read.delim(f, colClasses=c("character", "factor", "character",
                         "numeric", "numeric", "numeric", "numeric", "numeric",
                         "numeric", "numeric", "NULL"), strip.white=TRUE) } )
    names(DFs) <- paste("file", seq(along=DFs), sep="")  # name them

    ##  build shared set of variable names
    varNames <- Reduce(intersect, lapply(DFs, names))

    ## make sure that the data frames all have the correct variables
    wants <- c("Group", "Distance", "Aim.X", "Aim.Y", "Point.X", "Point.Y")
    has   <- wants %in% varNames
    if(!all(has)) {
        warning(c("At least one file is missing variable(s)\n",
                  paste(wants[!has], collapse= " "),
                  "\nthat may be required later by analysis functions"))
    }

    if(combine) { return(combineData(DFs)) } else { return(DFs) }
}

## read files from OnTarget-output version 2.*, 3.7*, 3.8*
readDataOT2 <-
function(fPath=".", fNames, fPat, combine=TRUE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1) { stop("No files were selected") }

    ## assumed variables: Project Title, Group, Ammunition, Distance,
    ## Aim X, Aim Y, Center X, Center Y, Point X, Point Y, Velocity (optional)
    ## 10 or 11 fields
    nFields <- unlist(lapply(files, function(x) count.fields(x, sep=",")))
    colClasses <- if(all(nFields == 10)) {
        c("character", "factor", "character", "numeric", "numeric",
          "numeric", "numeric", "numeric", "numeric", "numeric")
    } else if(all(nFields == 11)) {
        c("character", "factor", "character", "numeric", "numeric",
          "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
    } else {
        stop(c("It appears at least one file does not contain exactly\n",
               "the required set of variables - see help(readDataOT2)\n",
               "maybe you should use readDataMisc() instead"))
    }

    ## read in files into a list of data frames
    DFs <- lapply(files, function(f) {
        read.csv(f, colClasses=colClasses, strip.white=TRUE) } )
    names(DFs) <- paste("file", seq(along=DFs), sep="")  # name them

    ##  build shared set of variable names
    varNames <- Reduce(intersect, lapply(DFs, names))

    ## make sure that the data frames all have the correct variables
    wants <- c("Group", "Distance", "Aim.X", "Aim.Y", "Point.X", "Point.Y")
    has   <- wants %in% varNames
    if(!all(has)) {
        warning(c("At least one file is missing variable(s)\n",
                  paste(wants[!has], collapse= " "),
                  "\nthat may be required later by analysis functions"))
    }

    if(combine) {
        return(combineData(DFs))
    } else {
        return(DFs)
    }
}
