## convert size to angular measure
## MOA (= arcmin)
## Shooter's MOA (SMOA = IPHY, inches per hundred yards)
## milliradian milrad (1/1000 of a rad)
getMOA <-
function(x, dst, conversion="m2cm", type=c("MOA", "SMOA", "milrad")) {
    if(!is.numeric(x))   { stop(  "x must be numeric") }
    if(!is.numeric(dst)) { stop("dst must be numeric") }

    keep <- !is.na(x)
    if(any(x[keep] < 0)) { stop(  "x must be positive throughout") }
    if(any(dst <= 0))    { stop("dst must be positive throughout") }
    type <- match.arg(type)

    ## convert distance measure to the unit of the (x,y)-coordinates
    dstCommon <- getConvFac(conversion) * dst

    ## MOA
    ## calculate angle in rad, convert to deg, convert to MOA
    ## one degree = 1/360 of a circle's arc
    ## one degree = 60 MOA -> 1 MOA = 1/60 degree
    ## rad <- 2*atan(x/(2*dstCommon))  # size as angle in rad
    ## deg <- (360/(2*pi)) * rad       # size as angle in degree
    ## MOA <- 60*deg                   # degree in MOA

    ## SMOA
    ## get SMOA for 1 MOA = conversion factor for MOA2SMOA
    ## MOA2SMOA <- (pi/21600) * (1/atan(1/7200))
    MOA2SMOA       <- 1.0471975579301210602603044918909607385383783485554
    cnst21600pi    <- 6875.4935415698785052157785776926204398886566959877
    cnst1atan17200 <- 7200.0000462962960581466263965930916898334026996030

    ## milrad
    ## calculate angle in rad, convert to milliradian
    ## 1 rad = 1000 mils -> 1 milrad = 1/1000 rad
    ## rad <- atan(x/dstCommon)      # size as angle in rad
    ## milrad <- 1000*rad

    ## angle <- switch(type,
    ##     MOA =(21600/pi)*atan(x/(2*dstCommon)),  # size in arcmin
    ##     SMOA=(1/atan(1/7200))*atan(x/(2*dstCommon)),  # size in SMOA
    ##   milrad=2000*atan(x/(2*dstCommon)))        # size in milliradian
    angle <- x
    angle[keep] <- switch(type,
            MOA =cnst21600pi   *atan(x[keep]/(2*dstCommon)), # size in arcmin
            SMOA=cnst1atan17200*atan(x[keep]/(2*dstCommon)), # size in SMOA
          milrad=2000          *atan(x[keep]/(2*dstCommon))) # size in milliradian

    return(angle)
}

## convert angular measure to size
## MOA (= arcmin)
## Shooter's MOA (SMOA = IPHY, inches per hundred yards)
## milliradian milrad (1/1000 of a rad)
fromMOA <-
function(x, dst, conversion="m2cm", type=c("MOA", "SMOA", "milrad")) {
    if(!is.numeric(x))   { stop(  "x must be numeric") }
    if(!is.numeric(dst)) { stop("dst must be numeric") }
    keep <- !is.na(x)
    if(any(x[keep] < 0)) { stop(  "x must be positive throughout") }
    if(any(dst <= 0))    { stop("dst must be positive throughout") }
    type <- match.arg(type)

    ## convert distance measure to the unit of the (x,y)-coordinates
    dstCommon <- getConvFac(conversion) * dst

    ## MOA
    ## convert MOA to degree, convert to rad, calculate size
    ## one degree = 1/360 of a circle's arc
    ## one degree = 60 MOA -> 1 MOA = 1/60 degree
    ## deg <- MOA / 60                   # MOA as angle in degree
    ## rad <- (2*pi/360) * deg           # deg as angle in rad
    ## x   <- 2*dstCommon * tan(rad/2)   # angle in rad as size

    ## SMOA
    ## get SMOA for 1 MOA = conversion factor for MOA2SMOA
    ## MOA2SMOA <- (21600/pi)*atan(1/7200)
    MOA2SMOA    <- 0.95492965241113508367872462693595042591621952511092
    cnstPi21600 <- 0.00014544410433286079807697423070738439278690599071181

    ## milrad
    ## convert milliradian to rad, calculate size
    ## rad <- (1/1000)*milrad            # milliradian as rad
    ## x   <- 2*dstCommon * tan(rad/2)   # angle in rad as size

    ## size <- switch(type,
    ##     MOA =2*dstCommon*tan(x*pi/21600), # size from arcmin
    ##     SMOA=MOA2SMOA * 2*dstCommon*tan(x*pi/21600), # size from SMOA
    ##   milrad=2*dstCommon*tan(0.0005*x))   # size from milliradian
    size <- x
    size[keep] <- switch(type,
           MOA =2*dstCommon*tan(x[keep]*cnstPi21600), # size from arcmin
           SMOA=MOA2SMOA*2*dstCommon*tan(x[keep]*cnstPi21600),  # size from SMOA
         milrad=2*dstCommon*tan(0.0005*x[keep]))      # size from milliradian

    return(size)
}

makeMOA <-
function(x, dst, conversion) {
    if(length(x) == 1) {
        c(unit=x,
           MOA=getMOA(x, dst=dst, conversion=conversion, type="MOA"),
          SMOA=getMOA(x, dst=dst, conversion=conversion, type="SMOA"),
        milrad=getMOA(x, dst=dst, conversion=conversion, type="milrad"))
    } else {
        rbind(unit=x,
               MOA=getMOA(x, dst=dst, conversion=conversion, type="MOA"),
              SMOA=getMOA(x, dst=dst, conversion=conversion, type="SMOA"),
            milrad=getMOA(x, dst=dst, conversion=conversion, type="milrad"))
    }
}

## TODO: range estimation -> input size, MOA/SMOA/milrad, output distance
