getConvFac <-
function(conversion="m2cm") {
    ## check how the conversion factor is indicated
    convFac <- if(is.character(conversion)) {
        idxMM2MM <- conversion %in% c("mm2mm")
        idxMM2CM <- conversion %in% c("mm2cm")
        idxMM2M  <- conversion %in% c("mm2m")
        idxMM2YD <- conversion %in% c("mm2yd", "mm2yard", "mm2yards")
        idxMM2FT <- conversion %in% c("mm2ft", "mm2feet")
        idxMM2IN <- conversion %in% c("mm2in", "mm2inch")

        idxCM2MM <- conversion %in% c("cm2mm")
        idxCM2CM <- conversion %in% c("cm2cm")
        idxCM2M  <- conversion %in% c("cm2m")
        idxCM2YD <- conversion %in% c("cm2yd", "cm2yard", "cm2yards")
        idxCM2FT <- conversion %in% c("cm2ft", "cm2feet")
        idxCM2IN <- conversion %in% c("cm2in", "cm2inch", "cm2inches")

        idxM2MM  <- conversion %in% c("m2mm")
        idxM2CM  <- conversion %in% c("m2cm")
        idxM2M   <- conversion %in% c("m2m")
        idxM2YD  <- conversion %in% c("m2yd", "m2yard", "m2yards")
        idxM2FT  <- conversion %in% c("m2ft", "m2feet")
        idxM2IN  <- conversion %in% c("m2in", "m2inch", "m2inches")

        idxYD2MM <- conversion %in% c("yd2mm", "yard2mm", "yards2mm")
        idxYD2CM <- conversion %in% c("yd2cm", "yard2cm", "yards2cm")
        idxYD2M  <- conversion %in% c("yd2m",  "yard2m",  "yards2m")
        idxYD2YD <- conversion %in% c("yd2yd", "yard2yard", "yards2yards")
        idxYD2FT <- conversion %in% c("yd2ft", "yd2feet", "yard2ft", "yard2feet", "yards2ft", "yards2feet")
        idxYD2IN <- conversion %in% c("yd2in", "yard2inch", "yards2inch", "yards2inches")

        idxFT2MM <- conversion %in% c("ft2mm", "feet2mm")
        idxFT2CM <- conversion %in% c("ft2cm", "feet2cm")
        idxFT2M  <- conversion %in% c("ft2m",  "feet2m")
        idxFT2YD <- conversion %in% c("ft2yd", "feet2yd", "ft2yard", "feet2yard", "ft2yards", "feet2yards")
        idxFT2FT <- conversion %in% c("ft2ft", "feet2feet")
        idxFT2IN <- conversion %in% c("ft2in", "ft2inch", "ft2inches", "feet2in", "feet2inch", "feet2inches")
        
       ## did we catch all requested conversion types?
        idxAll <- idxMM2MM | idxMM2CM | idxMM2M | idxMM2YD | idxMM2FT | idxMM2IN |
                  idxCM2MM | idxCM2CM | idxCM2M | idxCM2YD | idxCM2FT | idxCM2IN |
                  idxM2MM  | idxM2CM  | idxM2M  | idxM2YD  | idxM2FT  | idxM2IN  |
                  idxYD2MM | idxYD2CM | idxYD2M | idxYD2YD | idxYD2FT | idxYD2IN |
                  idxFT2MM | idxFT2CM | idxFT2M | idxFT2YD | idxFT2FT | idxFT2IN
        if(!all(idxAll)) {
            warning(c('Conversion type(s) "', paste(conversion[!idxAll], collapse=", "),
                      '" not found - conversion factor set to 1'))
        }

        res <- rep(1, length(conversion))
        res[idxMM2MM] <-    1
        res[idxMM2CM] <-    0.1
        res[idxMM2M]  <-    0.001
        res[idxMM2YD] <-    0.0010936133
        res[idxMM2FT] <-    0.0032808399
        res[idxMM2IN] <-    0.0393700787

        res[idxCM2MM] <-   10
        res[idxCM2CM] <-    1
        res[idxCM2M]  <-    0.01
        res[idxCM2YD] <-    0.010936133
        res[idxCM2FT] <-    0.032808399
        res[idxCM2IN] <-    0.393700787

        res[ idxM2MM] <- 1000
        res[ idxM2CM] <-  100
        res[ idxM2M]  <-    1
        res[ idxM2YD] <-    1.0936133
        res[ idxM2FT] <-    3.2808399
        res[ idxM2IN] <-   39.3700787
        
        res[idxYD2MM] <-  914.4
        res[idxYD2CM] <-   91.44
        res[idxYD2M]  <-    0.9144
        res[idxYD2YD] <-    1
        res[idxYD2FT] <-    3
        res[idxYD2IN] <-   36

        res[idxFT2MM] <-  300.48
        res[idxFT2CM] <-   30.48
        res[idxFT2M]  <-    0.3048
        res[idxFT2YD] <-    0.333333333
        res[idxFT2FT] <-    1
        res[idxFT2IN] <-   12
        res
    } else if(is.numeric(conversion)) {  # factor is given directly
        if(any(conversion <= 0)) {
            stop("Conversion factor must be positive")
        }
        conversion
    } else {
        stop("Conversion must be a character constant or numeric")
    }                                # if(is.character(conversion))

    return(convFac)
}

# determine unit of (x,y)-coords from conversion string
getUnit <-
function(x="m2cm") {
    convEnd <- sub("^.+2(cm|mm|in|inch|inches)$", "\\1", x, ignore.case=TRUE, perl=TRUE)
    unit    <- substring(convEnd, 1, 2)
    if(!(unit %in% c("cm", "mm", "in"))) {
        stop("Unit not recognized - needs to have form *2cm, *2mm, *2in, *2inch, *2inches")
    } else {
        unit
    }
}

getMOA <-
function(x, dst, conversion="m2cm") {
    if(!is.numeric(x))   { stop(  "x must be numeric") }
    if(!is.numeric(dst)) { stop("dst must be numeric") }
    if(any(x < 0))       { stop(  "x must be positive throughout") }
    if(any(dst <= 0))    { stop("dst must be positive throughout") }

    ## get the conversion factor
    convFac <- getConvFac(conversion)

    ## convert distance measure to the unit of the (x,y)-coordinates
    dstCommon <- convFac * dst

    ## calculate angle
    # arc <- atan(x / dstCommon)         # size as angle in arc
    # deg <- arc*180/pi                  # size as angle in degree

    ## convert angle to MOA
    ## one degree = 1/360 of a circle's arc
    ## one degree = 60 MOA -> 1 MOA = 1/60 degree
    return(60*180*2 * atan(x / (2*dstCommon)) / pi)  # size in minutes of angle
}

fromMOA <-
function(x, dst, conversion="m2cm") {
    if(!is.numeric(x))   { stop(  "x must be numeric") }
    if(!is.numeric(dst)) { stop("dst must be numeric") }
    if(any(x < 0))       { stop(  "x must be positive throughout") }
    if(any(dst <= 0))    { stop("dst must be positive throughout") }

    ## get the conversion factor
    convFac <- getConvFac(conversion)

    ## convert distance measure to the unit of the (x,y)-coordinates
    dstCommon <- convFac * dst

    ## convert MOA to angle
    ## one degree = 1/360 of a circle's arc
    ## one degree = 60 MOA -> 1 MOA = 1/60 degree
    # deg <- x / 60                      # size as angle in degree
    # arc <- pi*deg/180                  # size as angle in arc

    return(tan(pi*x / (2*180*60)) * 2*dstCommon)  # size
}
