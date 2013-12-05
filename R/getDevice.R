## multiple calls to dev.new() currently don't work under RStudio
## do nothing instead
getDevice <-
function() {
    ## if under RStudio -> change device to explicit call
    dev <- getOption("device")
    if(class(dev) == "character") {
        if(dev == "RStudioGD") {
            dev <- function(...) {}
#             ## find out which function to use for opening a new plot window
#             ## find out operating system
#             osType <- .Platform$OS.type
#             if(osType == "windows") {
#                 dev <- windows
#             } else if(osType == "unix") {
#                 sysName <- Sys.info()[["sysname"]]
#                 if(sysName == "Linux") {
#                     dev <- x11
#                 } else if(sysName == "Darwin") {
#                     dev <- quartz
#                 }
#             }
        }
    }

    return(dev)
}
