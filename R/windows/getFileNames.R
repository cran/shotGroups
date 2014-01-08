getFileNames <-
function(fPath=getwd(), fNames, fPat) {
    files <- character(0)
    ## do we have file names or a name pattern?
    if(!missing(fNames)) {               # we have file names
        files <- paste(fPath, fNames, sep="/")
    } else if(!missing(fPat)) {          # we have a name pattern
        files <- list.files(path=fPath, pattern=fPat, full.names=TRUE)
    ## are we are in interactive mode AND under Windows?
	## we are under Windows since this sits in a platform-specific directory
#    } else if(interactive() && (.Platform$OS.type == "windows")) {
    } else if(interactive()) {
        ## choose files interactively
        myFilt <- rbind(Filters, txtCsvDat=c("Data files (*.txt, *.csv, *.dat)",
                                             "*.txt;*.csv;*.dat"))
        files <- choose.files(filters=myFilt[c("txtCsvDat", "All"), ], index=1)
    }

    return(files)
}
