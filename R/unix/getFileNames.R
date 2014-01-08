## under Unix-like OS, we don't have choose.files() and Filters
getFileNames <-
function(fPath=getwd(), fNames, fPat) {
    files <- character(0)
    ## do we have file names or a name pattern?
    if(!missing(fNames)) {               # we have file names
        files <- paste(fPath, fNames, sep="/")
    } else if(!missing(fPat)) {          # we have a name pattern
        files <- list.files(path=fPath, pattern=fPat, full.names=TRUE)
    }

    return(files)
}
