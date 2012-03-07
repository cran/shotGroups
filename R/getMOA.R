getConvFac <-
function(conversion="m2cm") {
    ## check how the conversion factor is indicated
    convFac <- if(is.character(conversion)) {
       idxM2CM  <- conversion %in% c("m2cm")
       idxM2MM  <- conversion %in% c("m2mm")
       idxYD2IN <- conversion %in% c("yd2in", "yard2inch", "yards2inch", "yards2inches")
       idxFT2IN <- conversion %in% c("ft2in", "ft2inch", "ft2inches", "feet2in",
                                     "feet2inch", "feet2inches")

       ## did we catch all conversion types?
       if(!all(idxM2CM | idxM2MM | idxYD2IN | idxFT2IN)) {
           warning("conversion type not found, conversion factor set to 1")
       }

       res <- numeric(length(conversion)) + 1
       res[idxM2CM]  <- 100
       res[idxM2MM]  <- 1000
       res[idxYD2IN] <-  36
       res[idxFT2IN] <-  12
       res
    } else if(is.numeric(conversion)) {  # factor is given directly
        if(any(conversion <= 0)) { stop("conversion factor must be positive") }
        conversion
    } else { stop("conversion must be a character constant or numeric") }

    return(convFac)
}

getMOA <-
function(x, dst, conversion="m2cm") {
    if(!is.numeric(x))   { stop("x must be numeric") }
    if(!is.numeric(dst)) { stop("dst must be numeric") }
    if(any(x < 0))       { stop("x must be positive throughout") }
    if(any(dst <= 0))    { stop("distance must be positive throughout") }

    ## get the conversion factor
    convFac <- getConvFac(conversion)

    ## convert distance measure to the unit of the (x,y)-coordinates
    dstCommon <- convFac * dst

    ## calculate angle
    # arc <- 2 * atan(x / (2*dstCommon))   # size as angle in arc
    # deg <- arc*180/pi                    # size as angle in degree

    ## convert angle to MOA
    ## one degree = 1/360 of a circle's arc
    ## one degree = 60 MOA -> 1 MOA = 1/60 degree
    return(60*360 * atan(x / (2*dstCommon)) / pi)  # size in minutes of angle
}

fromMOA <-
function(x, dst, conversion="m2cm") {
    if(!is.numeric(x))   { stop("x must be numeric") }
    if(!is.numeric(dst)) { stop("dst must be numeric") }
    if(any(x < 0))       { stop("x must be positive throughout") }
    if(any(dst <= 0))    { stop("distance must be positive throughout") }

    ## get the conversion factor
    convFac <- getConvFac(conversion)

    ## convert distance measure to the unit of the (x,y)-coordinates
    dstCommon <- convFac * dst

    ## convert MOA to angle
    ## one degree = 1/360 of a circle's arc
    ## one degree = 60 MOA -> 1 MOA = 1/60 degree
    # deg <- x / 60                        # size as angle in degree
    # arc <- pi*deg/180                    # size as angle in arc

    return(tan(pi*x / (180*60)) * dstCommon)  # size
}
