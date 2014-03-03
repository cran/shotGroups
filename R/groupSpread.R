groupSpread <-
function(xy, plots=TRUE, level=0.95, CEPtype="CorrNormal",
         bootCI=c("basic", "bca"), dstTarget=100, conversion="m2cm") {
    UseMethod("groupSpread")
}

groupSpread.data.frame <-
function(xy, plots=TRUE, level=0.95, CEPtype="CorrNormal",
         bootCI=c("basic", "bca"), dstTarget=100, conversion="m2cm") {
    xy <- getXYmat(xy)
    NextMethod("groupSpread")
}

groupSpread.default <-
function(xy, plots=TRUE, level=0.95, CEPtype="CorrNormal",
         bootCI=c("basic", "bca"), dstTarget=100, conversion="m2cm") {
    if(!is.matrix(xy))     { stop("xy must be a matrix") }
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(ncol(xy) != 2)      { stop("xy must have two columns") }
    if(!is.numeric(level)) { stop("level must be numeric") }
    if(level <= 0)         { stop("level must be > 0") }

    bootCI <- match.arg(bootCI, choices=c("none", "norm", "basic", "perc", "bca"), several.ok=TRUE)

    ## check if CI level is given in percent
    if(level >= 1) {
        while(level >= 1) { level <- level / 100 }
        warning(c("level must be in (0,1) and was set to ", level))
    }

    #####-----------------------------------------------------------------------
    ## prepare data
    X    <- xy[ , 1]                     # x-coords
    Y    <- xy[ , 2]                     # y-coords
    Npts <- nrow(xy)                     # number of observations
    res  <- vector("list", 0)            # empty list to later collect the results

    haveRob <- TRUE                      # can we do robust estimation?
    if(Npts < 4) {
        warning("We need >= 4 points for robust estimations")
        haveRob <- FALSE
    }                                    # if(haveRob)

    ## to determine axis limits later, collect all results in a vector
    axesCollX <- numeric(0)
    axesCollY <- numeric(0)

    #####-----------------------------------------------------------------------
    ## nonparametric bootstrap-CIs (basic and BCa)
    ## for all parameters where a CI is later reported
    if(!("none" %in% bootCI)) {          # do bootstrap CIs
        NrplMin <- 1499                  # minimum number of replications
        Nrpl <- if("bca" %in% bootCI) {  # number of replications
            max(NrplMin, Npts+1)         # BCa needs at least number of points
        } else {
            NrplMin
        }

        ## sdX, sdY, sigma, RSD, MR for one replication
        getSdSigRSDMR <- function(x, idx) {
            sdXY <- sqrt(diag(cov(x[idx, ])))
            rayParam <- getRayParam(x[idx, ], level=level, accuracy=FALSE)
            sigmaHat <- rayParam$sigma["sigma"]
            RSDhat   <- rayParam$RSD["RSD"]
            MRhat    <- rayParam$MR["MR"]
            return(c(sdXY, sigmaHat, RSDhat, MRhat))
        }

        bs <- boot::boot(xy, statistic=getSdSigRSDMR, R=Nrpl) # bootstrap

        ## extract CIs for all returned parameters
        sdXciBoot <- boot::boot.ci(bs, conf=level, type=bootCI, index=1)
        sdYciBoot <- boot::boot.ci(bs, conf=level, type=bootCI, index=2)
        sigCIboot <- boot::boot.ci(bs, conf=level, type=bootCI, index=3)
        RSDciBoot <- boot::boot.ci(bs, conf=level, type=bootCI, index=4)
         MRciBoot <- boot::boot.ci(bs, conf=level, type=bootCI, index=5)

        ## CI type names in output structure of boot.ci()
        CInames  <- c(basic="basic", norm="normal", perc="percent", bca="bca")
        CItype   <- CInames[bootCI]
        sdXciBoot <- c(sapply(CItype, function(x) {
            len  <- length(sdXciBoot[[x]])
            sdXciBoot[[x]][(len-1):len] }))
        sdYciBoot <- c(sapply(CItype, function(x) {
            len  <- length(sdYciBoot[[x]])
            sdYciBoot[[x]][(len-1):len] }))
        sigCIboot <- c(sapply(CItype, function(x) {
            len  <- length(sigCIboot[[x]])
            sigCIboot[[x]][(len-1):len] }))
        RSDciBoot <- c(sapply(CItype, function(x) {
            len  <- length(RSDciBoot[[x]])
            RSDciBoot[[x]][(len-1):len] }))
         MRciBoot <- c(sapply(CItype, function(x) {
            len  <- length(MRciBoot[[x]])
            MRciBoot[[x]][(len-1):len] }))

        names(sdXciBoot) <- paste("sdX",   rep(bootCI, each=2), c("(", ")"))
        names(sdYciBoot) <- paste("sdY",   rep(bootCI, each=2), c("(", ")"))
        names(sigCIboot) <- paste("sigma", rep(bootCI, each=2), c("(", ")"))
        names(RSDciBoot) <- paste("RSD",   rep(bootCI, each=2), c("(", ")"))
        names(MRciBoot)  <- paste("MR",    rep(bootCI, each=2), c("(", ")"))
    } else {
        sdXciBoot <- numeric(0)
        sdYciBoot <- numeric(0)
        sigCIboot <- numeric(0)
        RSDciBoot <- numeric(0)
         MRciBoot <- numeric(0)
    }

    #####-----------------------------------------------------------------------
    ## standard deviations of x- and y-coords
    varXY    <- diag(cov(xy))
    sdXY     <- sqrt(varXY)
    res$sdXY <- makeMOA(sdXY, dst=dstTarget, conversion=conversion)

    ## parametric CIs for true sd
    alpha <- 1-level
    sdXci <- sqrt((Npts-1)*varXY[1] / qchisq(c(1-(alpha/2), alpha/2), Npts-1))
    sdYci <- sqrt((Npts-1)*varXY[2] / qchisq(c(1-(alpha/2), alpha/2), Npts-1))
    names(sdXci) <- c("sdX (", "sdX )")
    names(sdYci) <- c("sdY (", "sdY )")

    ## combine parametric and bootstrap CIs and add to results
    res$sdXci <- makeMOA(c(sdXci, sdXciBoot), dst=dstTarget, conversion=conversion)
    res$sdYci <- makeMOA(c(sdYci, sdYciBoot), dst=dstTarget, conversion=conversion)

    ## robust standard deviations of x- and y-coords
    if(haveRob) {
        rob <- robustbase::covMcd(xy, cor=FALSE)
        sdXYrob <- sqrt(diag(rob$cov))
        res$sdXYrob <- makeMOA(sdXYrob, dst=dstTarget, conversion=conversion)
    } else {
        res$sdXYrob <- NULL
    }                                    # if(haveRob)

    ## (robust) center and covariance-matrix
    ctr <- colMeans(xy)                  # group center
    if(haveRob) {
        ctrRob <- rob$center
    }                                    # if(haveRob)
    res$covXY <- cov(xy)                 # covariance matrix
    if(haveRob) {
        res$covXYrob <- rob$cov
    } else {
        res$covXYrob <- NULL
    }                                    # if(haveRob)

    #####-----------------------------------------------------------------------
    ## mean distance to group center and associated parameterss
    dstCtr     <- getDistToCtr(xy)
    dstCtrMean <- mean(dstCtr)
    dstCtrMed  <- median(dstCtr)

    ## radial standard deviation
    ## http://ballistipedia.com/index.php?title=Describing_Precision
    rayParam <- getRayParam(xy, level=level, accuracy=FALSE)

    ## sigma, RSD, MR estimates with parametric confidence intervals
    sigma <- rayParam$sigma
    RSD   <- rayParam$RSD
    MR    <- rayParam$MR

    names(sigma) <- c("sigma", "sigma (", "sigma )")
    names(RSD)   <- c("RSD",   "RSD (",   "RSD )")
    names(MR)    <- c("MR",    "MR (",    "MR )")

    mDTCsigRSDmr  <- c(mean=dstCtrMean, median=dstCtrMed,
                       sigma["sigma"], RSD["RSD"], MR["MR"])
    res$distToCtr <- makeMOA(mDTCsigRSDmr, dst=dstTarget, conversion=conversion)

    ## combine parametric and bootstrap CIs and add to results
    sigmaCI     <- c(sigma[c("sigma (", "sigma )")], sigCIboot)
    res$sigmaCI <- makeMOA(sigmaCI, dst=dstTarget, conversion=conversion)

    RSDci     <- c(RSD[c("RSD (", "RSD )")], RSDciBoot)
    res$RSDci <- makeMOA(RSDci, dst=dstTarget, conversion=conversion)

    MRci     <- c(MR[c("MR (", "MR )")], MRciBoot)
    res$MRci <- makeMOA(MRci, dst=dstTarget, conversion=conversion)

    #####-----------------------------------------------------------------------
    ## maximum pairwise distance
    maxPD <- getMaxPairDist(xy)
    res$maxPairDist <- makeMOA(maxPD$d, dst=dstTarget, conversion=conversion)

    #####-----------------------------------------------------------------------
    ## width, height, FoM, diagonal of (minimum) bounding box
    bb    <- getBoundingBox(xy)          # bounding box
    bbMin <- getMinBBox(xy)              # minimum bounding box
    groupRect    <- c(width=bb$width,    height=bb$height,    FoM=bb$FoM,    diag=bb$diag)
    groupRectMin <- c(width=bbMin$width, height=bbMin$height, FoM=bbMin$FoM, diag=bbMin$diag)

    res$groupRect    <- makeMOA(groupRect,    dst=dstTarget, conversion=conversion)
    res$groupRectMin <- makeMOA(groupRectMin, dst=dstTarget, conversion=conversion)

    ## for axis limits
    axesCollX <- c(axesCollX, bb$pts[c(1, 3)], bbMin$pts[ , 1])
    axesCollY <- c(axesCollY, bb$pts[c(2, 4)], bbMin$pts[ , 2])

    #####-----------------------------------------------------------------------
    ## radius of minimum enclosing circle
    minCirc <- getMinCircle(xy)          # minimum enclosing circle
    res$minCircleRad <- makeMOA(minCirc$rad, dst=dstTarget, conversion=conversion)

    ## for axis limits
    axesCollX <- c(axesCollX, minCirc$ctr[1] + minCirc$rad,
                              minCirc$ctr[1] - minCirc$rad)
    axesCollY <- c(axesCollY, minCirc$ctr[2] + minCirc$rad,
                              minCirc$ctr[2] - minCirc$rad)

    #####-----------------------------------------------------------------------
    ## confidence ellipse measures
    confEll <- getConfEll(xy, level, dstTarget, conversion, doRob=haveRob)
    res$confEll <- confEll$size

    ## for axis limits
    axesCollX <- c(axesCollX, confEll$ctr[1] + confEll$size["unit", "semi-major"],
                              confEll$ctr[1] - confEll$size["unit", "semi-major"])
    axesCollY <- c(axesCollY, confEll$ctr[2] + confEll$size["unit", "semi-major"],
                              confEll$ctr[2] - confEll$size["unit", "semi-major"])

    if(haveRob) {
        res$confEllRob <- confEll$sizeRob

        ## for axis limits
        axesCollX <- c(axesCollX, confEll$ctrRob[1] + confEll$sizeRob["unit", "semi-major"],
                                  confEll$ctrRob[1] - confEll$sizeRob["unit", "semi-major"])
        axesCollY <- c(axesCollY, confEll$ctrRob[2] + confEll$sizeRob["unit", "semi-major"],
                                  confEll$ctrRob[2] - confEll$sizeRob["unit", "semi-major"])
    } else {
        res$confEllRob <- NULL
    }                                    # if(haveRob)

    res$confEllShape <- confEll$shape
    if(haveRob) {
        res$confEllShapeRob <- confEll$shapeRob
    } else {
        res$confEllShapeRob <- NULL
    }                                    # if(haveRob)

    #####-----------------------------------------------------------------------
    ## circular error probable
    CEP <- getCEP(xy, level=level, type=CEPtype, dstTarget=dstTarget, conversion=conversion, accuracy=FALSE)
    res$CEP <- CEP$CEP

    #####-----------------------------------------------------------------------
    ## plotting
    if(plots) {
        ## infer (x,y)-coord units from conversion
        unitXY  <- getUnits(conversion, first=FALSE)
        unitDst <- getUnits(conversion, first=TRUE)
        devNew  <- getDevice()           # platform-dependent window open

        #####-------------------------------------------------------------------
        ## diagram: histogram for distances to group center
        devNew()                         # open new diagram
        hist(dstCtr, breaks="FD", freq=FALSE,
             main="Histogram distances to center w/ kernel density estimate",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("distance [", unitXY, " ]"))
        rug(jitter(dstCtr))              # show single values

        ## add Rayleigh fit and kernel density estimate
        dRaySigma <- function(x) { dRayleigh(x, sigma["sigma"]) }
        curve(dRaySigma, lwd=2, col="blue", add=TRUE)
        lines(density(dstCtr), lwd=2, col="red")  # kernel density estimate
        legend(x="topright",
               legend=c("Rayleigh distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=1, lwd=2, bg=rgb(1, 1, 1, 0.7))

        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        ## determine axis limits
        xLims <- range(c(X, axesCollX))
        yLims <- range(c(Y, axesCollY))

        devNew()                           # open new diagram
        plot(Y ~ X, asp=1, xlim=xLims, ylim=yLims, pch=20,
             main="Group (x,y)-coordinates",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim

        ## add group center and robust estimate for group center
        points(ctr[1], ctr[2], col="red",  pch=4, lwd=2, cex=1.5)
        if(haveRob) {
            points(ctrRob[1], ctrRob[2], col="blue", pch=4, lwd=2, cex=1.5)
        }                                # if(haveRob)

        ## add confidence ellipses (parametric, robust),
        ## and a circle with mean distance to center
        drawEllipse(confEll, pch=4, fg="red", lwd=2)
        if(haveRob) {
            drawEllipse(ctrRob, res$covXYrob, radius=confEll$magFac,
                        pch=4, fg="blue", lwd=2)
        }                                # if(haveRob)
        drawCircle(ctr, radius=mean(dstCtr), fg="black", lwd=2)

        ## add legend
        legend(x="bottomleft", legend=c("center", "center (robust)",
               paste(100*level, "% confidence ellipse", sep=""),
               paste(100*level, "% confidence ellipse (robust)", sep=""),
               "mean distance to center"),
               col=c("red", "blue", "red", "blue", "black"),
               pch=c(4, 4, NA, NA, NA),
               lty=c(NA, NA, 1, 1, 1), lwd=2, bg=rgb(1, 1, 1, 0.7))

        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(Y ~ X, asp=1, xlim=xLims, ylim=yLims, pch=20,
             main="Group (x,y)-coordinates",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")                            # add point of aim
        points(ctr[1], ctr[2], col="gray40", pch=4, lwd=4, cex=2.5)  # add group center

        ## add bounding box, minimum bounding box, minimum enclosing circle,
        ## and maximum group spread
        drawBox(bb, fg="magenta", lwd=2)
        drawBox2(bbMin, fg="red", lwd=2)
        drawCircle(minCirc, fg="blue", lwd=2)
        segments(x0=xy[maxPD$idx[1], 1], y0=xy[maxPD$idx[1], 2],
                 x1=xy[maxPD$idx[2], 1], y1=xy[maxPD$idx[2], 2],
                 col="green3", lwd=2)

        ## add legend
        legend(x="bottomleft", legend=c("group center", "bounding box",
               "minimum bounding box",
               "minimum enclosing circle", "maximum group spread"),
               col=c("gray40", "magenta", "red", "blue", "green3"),
               pch=c(4, NA, NA, NA, NA), lty=c(NA, 1, 1, 1, 1), lwd=2,
               bg=rgb(1, 1, 1, 0.7))
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
