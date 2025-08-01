groupLocation <-
function(xy, level=0.95, plots=TRUE, bootCI="none",
         dstTarget, conversion) {
    UseMethod("groupLocation")
}

groupLocation.data.frame <-
function(xy, level=0.95, plots=TRUE, bootCI="none",
         dstTarget, conversion) {
    ## distance to target from override or from data
    if(missing(dstTarget)) {
        dstTarget <- if(hasName(xy, "distance")) {
            xy[["distance"]]
        } else {
            NA_real_
        }
    }
    
    ## determine conversion factor from data if override is not given
    if(missing(conversion)) {
        conversion <- determineConversion(xy)
    }

    xy <- getXYmat(xy)
    
    groupLocation(xy, level=level, plots=plots, bootCI=bootCI,
                  dstTarget=dstTarget, conversion=conversion)
}

groupLocation.default <-
function(xy, level=0.95, plots=TRUE, bootCI="none",
         dstTarget, conversion) {
    if(!is.matrix(xy))     { stop("xy must be a matrix") }
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(ncol(xy) != 2L)     { stop("xy must have two columns") }
    if(!is.numeric(level)) { stop("level must be numeric") }
    if(level <= 0)         { stop("level must be > 0") }

    bootCI <- match.arg(bootCI, choices=c("none", "norm", "basic", "perc", "bca"), several.ok=TRUE)

    ## check if CI level is given in percent
    if(level >= 1) {
        while(level >= 1) { level <- level / 100 }
        warning(c("level must be in (0,1) and was set to ", level))
    }

    dstTarget <- if(missing(dstTarget)    ||
                    all(is.na(dstTarget)) ||
                    (length(unique(dstTarget)) > 1L)) {
        NA_real_
    } else {
        mean(dstTarget)
    }
    
    conversion <- if(missing(conversion)    ||
                     all(is.na(conversion)) ||
                     (length(unique(conversion)) > 1L)) {
        NA_character_
    } else {
        unique(conversion)
    }
    
    #####-----------------------------------------------------------------------
    ## prepare data
    X    <- xy[ , 1]                     # x-coords
    Y    <- xy[ , 2]                     # y-coords
    Npts <- nrow(xy)                     # number of observations
    res  <- vector("list", 0)            # empty list to later collect the results

    ## can we do robust estimation?
    haveRobustbase <- requireNamespace("robustbase", quietly=TRUE)
    haveRob <- if(haveRobustbase && (Npts >= 4L)) {
        TRUE
    } else {
        if(Npts < 4L) {
            warning("We need >= 4 points for robust estimations")
        }
        
        if(!haveRobustbase) {
            warning("Please install package 'robustbase' for robust estimations")
        }
        
        FALSE
    }
    
    #####-----------------------------------------------------------------------
    ## location measures
    res$ctr <- colMeans(xy)              # center of joint (x,y)-distribution

    ## robust estimation of center
    res$ctrRob <- if(haveRob) {
        robustbase::covMcd(xy)$center
    } else {
        NULL
    }                                    # if(haveRob)

    distPOA     <- sqrt(sum(res$ctr^2))  # distance to point of aim
    res$distPOA <- makeMOA(distPOA, dst=dstTarget, conversion=conversion)

    res$distPOArob <- if(haveRob) {      # rob distance to point of aim
        distPOArob <- sqrt(sum(res$ctrRob^2))
        makeMOA(distPOArob, dst=dstTarget, conversion=conversion)
    } else {
        NULL
    }                                    # if(haveRob)

    ## Hotelling's T^2 test for equality of (x,y)-center with point of aim (0,0)
    res$Hotelling <- if(Npts > 2L) {
        anova(lm(cbind(X, Y) ~ 1), test="Hotelling-Lawley")
    } else {
        warning("We need >= 3 points for Hotelling's T^2 test")
        NULL
    }                                    # if(Npts > 2L)

    #####-----------------------------------------------------------------------
    ## confidence intervals for x- and y-coords
    ## parametric: t-CI
    alpha  <- 1-level                    # alpha-level
    tCrit  <- qt(c(alpha/2, 1-alpha/2), Npts-1)  # critical t-values left and right
    Mx     <- mean(X)                    # mean x-coords
    My     <- mean(Y)                    # mean y-coords
    sMx    <- sd(X) / sqrt(Npts)         # standard error of the mean x
    sMy    <- sd(Y) / sqrt(Npts)         # standard error of the mean y
    ctrXci <- rbind(t=rev(Mx-tCrit*sMx)) # t-CI x-coords
    ctrYci <- rbind(t=rev(My-tCrit*sMy)) # t-CI y-coords

    ## non-parametric: bootstrap-CIs for center (basic and BCa)
    if(!("none" %in% bootCI)) {          # do bootstrap CIs
        NrplMin <- 1499L                 # minimum number of replications
        Nrpl <- if("bca" %in% bootCI) {  # number of replications
            max(NrplMin, Npts+1)         # BCa needs at least this number of points
        } else {
            NrplMin
        }

        ## group center for one replication
        getCtr  <- function(x, idx) { colMeans(x[idx, , drop=FALSE]) }
        bs      <- boot::boot(xy, statistic=getCtr, R=Nrpl)  # bootstrap centers
        xCIboot <- boot::boot.ci(bs, conf=level, type=bootCI, index=1) # x
        yCIboot <- boot::boot.ci(bs, conf=level, type=bootCI, index=2) # y

        ## CI type names in output structure of boot.ci()
        CInames <- c(basic="basic", norm="normal", perc="percent", bca="bca")
        CItype  <- CInames[bootCI]
        xCImat  <- vapply(CItype, function(x) {
            len <- length(xCIboot[[x]])
            xCIboot[[x]][(len-1):len] }, numeric(2))
        yCImat  <- vapply(CItype, function(x) {
            len <- length(yCIboot[[x]])
            yCIboot[[x]][(len-1):len] }, numeric(2))

        ## add bootstrap CIs to parametric CI
        ctrXci <- rbind(ctrXci, t(xCImat))
        ctrYci <- rbind(ctrYci, t(yCImat))
    }

    res$ctrXci <- ctrXci
    res$ctrYci <- ctrYci
    colnames(res$ctrXci) <- c("x (", "x )")
    colnames(res$ctrYci) <- c("y (", "y )")

    if(plots) {
        ## infer (x,y)-coord units from conversion
        unitXY  <- na.omit(getUnits(conversion, first=FALSE))
        unitDst <- na.omit(getUnits(conversion, first=TRUE))
        devNew  <- getDevice()           # platform-dependent window open

        ## distance to target may be heterogeneous
        dstTargetPlot <- paste(unique(round(na.omit(dstTarget))), collapse=", ")
        
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(Y ~ X, asp=1, main="Group (x,y)-coordinates", pch=16,
             sub=paste("distance:", dstTargetPlot, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="gray")     # add point of aim

        ## add (robust) group center
        points(res$ctr[1], res$ctr[2], col="blue", pch=4, lwd=2, cex=1.5)

        if(haveRob) {
            points(res$ctrRob[1], res$ctrRob[2], col="red",
                   pch=4, lwd=2, cex=1.5)
        }                                # if(haveRob)

        ## add legend
        legend(x="bottomleft", legend=c("group center", "robust group center"),
               col=c("blue", "red"), pch=4, lty=NA, lwd=2, bg=rgb(1, 1, 1, 0.7))
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
