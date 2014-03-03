simRingCount <-
function(xy, target, caliber, unit="cm") {
    UseMethod("simRingCount")
}

simRingCount.data.frame <-
    function(xy, target, caliber, unit="cm") {
        xy <- getXYmat(xy)
    NextMethod("simRingCount")
}

simRingCount.default <-
function(xy, target, caliber, unit="cm") {
    if(!is.matrix(xy))       { stop("xy must be a matrix") }
    if(!is.numeric(xy))      { stop("xy must be numeric") }
    if(ncol(xy) != 2)        { stop("xy must have two columns") }
    if(!is.numeric(caliber)) { stop("caliber must be numeric") }
    if(caliber <= 0)         { stop("caliber must > 0") }

    unit <- match.arg(unit, choices=c("cm", "mm", "m", "in", "ft", "yd"))

    ## prepare data
    ## get target definition in requested unit
    trgt <- getTarget(target, unit=unit)

    ## convert caliber to required unit
    convFac <- getConvFac(paste0("mm2", unit))
    calSize <- convFac * caliber

    ## get distance of inner edge of bullet hole to point of aim (0,0)
    ## negative difference -> distance from other side
    dstPOA <- abs(sqrt(rowSums(xy^2)) - calSize/2)

    ## cut with breaks = ring radii
    rings <- if(!is.null(trgt$countMouche) && trgt$countMouche) {
        ## with 1st ring (mouche, ring 10 inner sub-division)
        with(trgt, cut(dstPOA, breaks=c(0, ringRu, Inf),
                       labels=c((maxVal+1):(maxVal-nRings+1), 0)),
                       include.lowest=TRUE)
    } else {
        ## except 1st ring (mouche, ring 10 inner sub-division)
        with(trgt, cut(dstPOA, breaks=c(0, ringRu[-1], Inf),
                       labels=c(maxVal:(maxVal-nRings+1), 0)),
                       include.lowest=TRUE)
    }

    ## convert factor labels to numeric, then sum
    ringCount <- sum(as.numeric(as.character(rings)))  # observed ring count
    ringMax   <- 10 * nrow(xy)                         # maximum possible

    return(list(count=ringCount, max=ringMax, rings=rings))
}
