groupSpread <-
function(xy, plots=TRUE, level=0.95,
         CEPtype=c("Rayleigh", "Grubbs", "RAND"),
         sigmaType=c("Rayleigh", "Gauss"),
         dstTarget=100, conversion="m2cm") {
    UseMethod("groupSpread")
}

groupSpread.data.frame <-
function(xy, plots=TRUE, level=0.95,
         CEPtype=c("Rayleigh", "Grubbs", "RAND"),
         sigmaType=c("Rayleigh", "Gauss"),
         dstTarget=100, conversion="m2cm") {
    xy <- getXYmat(xy)
    NextMethod("groupSpread")
}

groupSpread.default <-
function(xy, plots=TRUE, level=0.95,
         CEPtype=c("Rayleigh", "Grubbs", "RAND"),
         sigmaType=c("Rayleigh", "Gauss"),
         dstTarget=100, conversion="m2cm") {
    if(!is.matrix(xy))     { stop("xy must be a matrix") }
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(ncol(xy) != 2)      { stop("xy must have two columns") }
    if(!is.numeric(level)) { stop("level must be numeric") }
    if(level <= 0)         { stop("level must be > 0") }

    sigmaType <- match.arg(sigmaType)

    ## check if CI level is given in percent
    if(level > 1) {
        while(level > 1) {
            level <- level / 100
        }
        warning(c("level must be in (0,1) and was set to ", level, "\n"))
    }

    N <- nrow(xy)
    alpha <- 1-level

    haveRob <- TRUE                      # can we do robust estimation?
    if(nrow(xy) < 4) {
        warning("We need >= 4 points for robust estimations")
        haveRob <- FALSE
    }                                    # if(haveRob)

    #####-----------------------------------------------------------------------
    ## prepare data
    X   <- xy[ , 1]                      # x-coords
    Y   <- xy[ , 2]                      # y-coords
    res <- vector("list", 0)             # empty list to later collect the results

    #####-----------------------------------------------------------------------
    ## nonparametric bootstrap-CIs (basic and BCa)
    ## for all parameters where a CI is later reported
    nR <- 1499                           # number of replications
    ## sdX, sdY, sigma, RSD, MR for one replication
    getSdSigRSDMR <- function(x, idx) {
        sdXY <- sqrt(diag(cov(x[idx, ])))
        rayParam <- getRayParam(x[idx, ], type=sigmaType, level=level, accuracy=FALSE)
        sigmaHat <- rayParam$sigma["sigma"]
        RSDhat   <- rayParam$RSD["RSD"]
        MRhat    <- rayParam$MR["MR"]
        return(c(sdXY, sigmaHat, RSDhat, MRhat))
    }

    bs <- boot::boot(xy, statistic=getSdSigRSDMR, R=nR) # bootstrap

    ## extract CIs for all returned parameters
    sdXciBoot <- boot::boot.ci(bs, conf=level, type=c("basic", "bca"), index=1)
    sdYciBoot <- boot::boot.ci(bs, conf=level, type=c("basic", "bca"), index=2)
    sigCIboot <- boot::boot.ci(bs, conf=level, type=c("basic", "bca"), index=3)
    RSDciBoot <- boot::boot.ci(bs, conf=level, type=c("basic", "bca"), index=4)
     MRciBoot <- boot::boot.ci(bs, conf=level, type=c("basic", "bca"), index=5)
    
    sdXciBoot <- c(sdXciBoot$basic[4:5], sdXciBoot$bca[4:5])
    sdYciBoot <- c(sdYciBoot$basic[4:5], sdYciBoot$bca[4:5])
    sigCIboot <- c(sigCIboot$basic[4:5], sigCIboot$bca[4:5])
    RSDciBoot <- c(RSDciBoot$basic[4:5], RSDciBoot$bca[4:5])
    MRciBoot  <- c( MRciBoot$basic[4:5],  MRciBoot$bca[4:5])
    
    names(sdXciBoot) <- c(  "sdX basic (",   "sdX basic )",   "sdX BCa (",   "sdX BCa )")
    names(sdYciBoot) <- c(  "sdY basic (",   "sdY basic )",   "sdY BCa (",   "sdY BCa )")
    names(sigCIboot) <- c("sigma basic (", "sigma basic )", "sigma BCa (", "sigma BCa )")
    names(RSDciBoot) <- c(  "RSD basic (",   "RSD basic )",   "RSD BCa (",   "RSD BCa )")
    names(MRciBoot)  <- c(   "MR basic (",    "MR basic )",    "MR BCa (",    "MR BCa )")
    
    #####-----------------------------------------------------------------------
    ## standard deviations of x- and y-coords
    res$sdXY <- rbind(unit=c(sd(X), sd(Y)),      # standard deviations
                       MOA=getMOA(c(sd(X), sd(Y)), dst=dstTarget, conversion=conversion))
    colnames(res$sdXY) <- c("X", "Y")

    ## parametric CIs for true sd
    sdXci <- sqrt((N-1)*var(X) / qchisq(c(1-(alpha/2), alpha/2), N-1))
    sdYci <- sqrt((N-1)*var(Y) / qchisq(c(1-(alpha/2), alpha/2), N-1))
    names(sdXci) <- c("sdX (", "sdX )")
    names(sdYci) <- c("sdY (", "sdY )")

    ## combine parametric and bootstrap CIs and add to results
    res$sdXci <- rbind(unit=c(sdXci, sdXciBoot),
                        MOA=getMOA(c(sdXci, sdXciBoot), dst=dstTarget, conversion=conversion))
    res$sdYci <- rbind(unit=c(sdYci, sdYciBoot),
                        MOA=getMOA(c(sdYci, sdYciBoot), dst=dstTarget, conversion=conversion))

    ## robust standard deviations of x- and y-coords
    if(haveRob) {
        rob <- robustbase::covMcd(xy, cor=FALSE)
        res$sdXYrob <- rbind(unit=sqrt(diag(rob$cov)),   # robust sd
                              MOA=getMOA(sqrt(diag(rob$cov)), dst=dstTarget, conversion=conversion))
        colnames(res$sdXYrob) <- c("X", "Y")
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
    rayParam <- getRayParam(xy, type=sigmaType, level=level, accuracy=FALSE)

    ## sigma, RSD, MR estimates with parametric confidence intervals
    sigma <- rayParam$sigma
    RSD   <- rayParam$RSD
    MR    <- rayParam$MR

    names(sigma) <- c("sigma", "sigma (", "sigma )")
    names(RSD)   <- c("RSD",   "RSD (",   "RSD )")
    names(MR)    <- c("MR",    "MR (",    "MR )")

    mDTCsigRSDmr  <- c(mean=dstCtrMean, median=dstCtrMed,
                       sigma["sigma"], RSD["RSD"], MR["MR"])
    res$distToCtr <- rbind(unit=mDTCsigRSDmr,
                            MOA=getMOA(mDTCsigRSDmr, dst=dstTarget, conversion=conversion))

    ## combine parametric and bootstrap CIs and add to results
    res$sigmaCI <- rbind(unit=c(sigma[c("sigma (", "sigma )")], sigCIboot),
                          MOA=getMOA(c(sigma[c("sigma (", "sigma )")], sigCIboot), dst=dstTarget, conversion=conversion))
    
    res$RSDci <- rbind(unit=c(RSD[c("RSD (", "RSD )")], RSDciBoot),
                        MOA=getMOA(c(RSD[c("RSD (", "RSD )")], RSDciBoot), dst=dstTarget, conversion=conversion))

    res$MRci <- rbind(unit=c(MR[c("MR (", "MR )")], MRciBoot),
                       MOA=getMOA(c(MR[c("MR (", "MR )")], MRciBoot), dst=dstTarget, conversion=conversion))

    #####-----------------------------------------------------------------------
    ## maximum pairwise distance
    maxPD <- getMaxPairDist(xy)
    res$maxPairDist <- c(unit=maxPD$d,
                          MOA=getMOA(maxPD$d, dst=dstTarget, conversion=conversion))

    #####-----------------------------------------------------------------------
    ## width and height of (minimum) bounding box
    bb    <- getBoundingBox(xy)          # bounding box
    bbMin <- getMinBBox(xy)              # minimum bounding box
    groupWidth     <- c(bb$width,     getMOA(bb$width,     dst=dstTarget, conversion=conversion))
    groupHeight    <- c(bb$height,    getMOA(bb$height,    dst=dstTarget, conversion=conversion))
    groupWidthMin  <- c(bbMin$width,  getMOA(bbMin$width,  dst=dstTarget, conversion=conversion))
    groupHeightMin <- c(bbMin$height, getMOA(bbMin$height, dst=dstTarget, conversion=conversion))
    FoM            <- c(bb$FoM,       getMOA(bb$FoM,       dst=dstTarget, conversion=conversion))
    FoMmin         <- c(bbMin$FoM,    getMOA(bbMin$FoM,    dst=dstTarget, conversion=conversion))
    bbDiag         <- c(bb$diag,      getMOA(bb$diag,      dst=dstTarget, conversion=conversion))
    bbMinDiag      <- c(bbMin$diag,   getMOA(bbMin$diag,   dst=dstTarget, conversion=conversion))

    res$groupRect    <- cbind(width=groupWidth,    height=groupHeight,    FoM=FoM,    diag=bbDiag)
    res$groupRectMin <- cbind(width=groupWidthMin, height=groupHeightMin, FoM=FoMmin, diag=bbMinDiag)
    
    rownames(res$groupRect)    <- c("unit", "MOA")
    rownames(res$groupRectMin) <- c("unit", "MOA")

    #####-----------------------------------------------------------------------
    ## radius of minimum enclosing circle
    minCirc <- getMinCircle(xy)          # minimum enclosing circle
    res$minCircleRad <- c(unit=minCirc$rad,
                           MOA=getMOA(minCirc$rad, dst=dstTarget, conversion=conversion))

    #####-----------------------------------------------------------------------
    ## confidence ellipse measures
    confEll <- getConfEll(xy, level, dstTarget, conversion, doRob=haveRob)
    res$confEll <- confEll$size
    if(haveRob) {
        res$confEllRob <- confEll$sizeRob
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
    CEP <- getCEP(xy, dstTarget=dstTarget, conversion=conversion, accuracy=FALSE)
    res$CEP <- CEP[CEPtype]

    #####-----------------------------------------------------------------------
    ## plotting
    if(plots) {
        devNew <- getDevice()            # platform-dependent window open
        #####-------------------------------------------------------------------
        ## diagram: histogram for distances to group center
        devNew()                         # open new diagram
        hist(dstCtr, breaks="FD", freq=FALSE,
             main="Histogram distances to center w/ kernel density estimate")
        rug(jitter(dstCtr))              # show single values

        ## add Rayleigh fit and kernel density estimate
        dRaySigma <- function(x) { dRayleigh(x, sigma["sigma"]) }
        curve(dRaySigma, lwd=2, col="blue", add=TRUE)
        lines(density(dstCtr), lwd=2, col="red")  # kernel density estimate
        legend(x="topright",
               legend=c("Rayleigh distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=1, lwd=2, bg="white")
        
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        ## determine axis limits
        xLims <- range(c(bbMin$pts[ , 1], minCirc$ctr[1] + c(-minCirc$rad, minCirc$rad)))
        yLims <- range(c(bbMin$pts[ , 2], minCirc$ctr[2] + c(-minCirc$rad, minCirc$rad)))

        devNew()                           # open new diagram
        plot(Y ~ X, asp=1, xlim=xLims, ylim=yLims, pch=20,
             main="Group (x,y)-coordinates")
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
               lty=c(NA, NA, 1, 1, 1), lwd=2, bg="white")
        
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(Y ~ X, asp=1, xlim=xLims, ylim=yLims, pch=20,
             main="Group (x,y)-coordinates")
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
               pch=c(4, NA, NA, NA, NA), lty=c(NA, 1, 1, 1, 1), lwd=2, bg="white")
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
