## find out which function to use for opening a new plot window
## multiple calls to dev.new() currently don't work under RStudio
getDevice <-
function() {
    ## if under RStudio -> change device to explicit call
    dev <- getOption("device")
    if(class(dev) == "character") {
        if(dev == "RStudioGD") {                  # yes
            ## find out operating system
            osType <- .Platform$OS.type
            if(osType == "windows") {
                dev <- windows
            } else if(osType == "unix") {
                sysName <- Sys.info()[["sysname"]]
                if(sysName == "Linux") {
                    dev <- x11
                } else if(sysName == "Darwin") {
                    dev <- quartz
                }
            }
        }
    }

    return(dev)
}
