compareGroups <-
function(DF, plots=TRUE, xyTopLeft=TRUE, ABalt=c("two.sided", "less", "greater"),
         Walt=c("two.sided", "less", "greater"), CEPtype="CorrNormal", conversion="m2cm") {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    #####-----------------------------------------------------------------------
    ## make sure DF has the required variable names and at least two groups
    varNames <- names(DF)                # what variables are present
    needsSer <- "Series"                 # required
    needsXY1 <- c("Point.X", "Point.Y")  # coordinates must have this name
    needsXY2 <- c("X", "Y")              # or this
    wantsDst <- "Distance"               # useful
    wantsAIM <- c("Aim.X", "Aim.Y")      # useful
    hasSer   <- needsSer %in% varNames   # required ones we have
    hasXY1   <- needsXY1 %in% varNames   # coordinates we have
    hasXY2   <- needsXY2 %in% toupper(varNames)
    hasDst   <- wantsDst %in% varNames   # useful ones we have
    hasAIM   <- wantsAIM %in% varNames   # useful ones we have

    if(!all(hasSer)) {
        stop(c("The data frame is missing variable\n",
               paste(needsSer[!hasSer], collapse=" ")))
    }

    if(!xor(all(hasXY1), all(hasXY2))) {
        stop("Coordinates must be named either X, Y or Point.X, Point.Y")
    }

    if(!all(hasDst)) {
        warning(c("The data frame is missing variable\n",
                  paste(wantsDst[!hasDst], collapse=" "), "\n",
                  "Distance is assumed to be 100"))
        DF$Distance <- 100
    }

    if(!all(hasAIM)) {
        warning(c("The data frame is missing variable(s)\n",
                  paste(wantsAIM[!hasAIM], collapse=" "), "\n",
                  "Point of Aim is assumed to be in (0,0)"))
        DF$Aim.X <- 0
        DF$Aim.Y <- 0
    }

    #####-----------------------------------------------------------------------
    ## prepare data
    res <- vector("list", 0)             # empty list to later collect the results
    if(!is.factor(DF$Series)) {          # make sure Series is a factor
        DF$Series <- as.factor(DF$Series)
    } else {
        DF$Series <- droplevels(DF$Series) # remove all non-used factor levels
    }

    ## check if we have enough groups and points per group
    if(nlevels(DF$Series) < 2) { stop("We need >= 2 groups for a comparison") }
    if(any(xtabs(~ Series, data=DF) < 2)) { stop("We need >= 2 points in each group") }

    ## prepare data: get (x,y)-coords relative to point of aim as matrix
    xy <- getXYmat(DF, xyTopLeft=xyTopLeft)
    DF <- cbind(DF, xy)
    dstTarget <- tapply(DF$Distance, DF$Series, mean)  # distances to target

    ## for each group extract the new (x,y)-coords as a matrix
    extractXY <- function(x) {
        DFsub <- x[ , which(!(names(x) %in% c("X", "Y")))]
        getXYmat(DFsub, xyTopLeft=xyTopLeft)
    }

    xyL  <- lapply(split(DF, DF$Series), extractXY)
    nS   <- length(xyL)                  # total number of series
    nObs <- sapply(xyL, length)          # number of obs per series
    names(xyL) <- paste("Series", levels(DF$Series), sep="")

    ## to determine axis limits later, collect all results in a vector
    axisLimsX <- numeric(0)
    axisLimsY <- numeric(0)

    #####-----------------------------------------------------------------------
    ## location measures
    res$ctr <- sapply(xyL, colMeans)     # group centers
    distPOA <- sqrt(colSums(res$ctr^2))  # distances to point of aim
    distPOAmoa  <- Map(makeMOA, distPOA, dst=dstTarget, conversion=conversion)
    res$distPOA <- do.call("cbind", distPOAmoa)

    ## multivariate location test for equal group centers (relative to POA)
    res$MANOVA <- anova(lm(cbind(X, Y) ~ Series, data=DF), test="Wilks")

    #####-----------------------------------------------------------------------
    ## shape measures
    ## correlation matrices for x- and y-coords
    res$corXY <- lapply(xyL, cor)

    #####-----------------------------------------------------------------------
    ## spread measures
    ## standard deviations for x- and y-coords with parametric CIs
    alpha <- 0.05
    getSDxyCI <- function(x) {
        N    <- nrow(x)
        sdXY <- sqrt(diag(cov(x)))       # standard deviations

        ## and their parametric 95%-confidence intervals
        sdXci  <- sqrt((N-1)*sdXY["X"]^2 / qchisq(c(1-(alpha/2), alpha/2), N-1))
        sdYci  <- sqrt((N-1)*sdXY["Y"]^2 / qchisq(c(1-(alpha/2), alpha/2), N-1))
        sdXYci <- c(sdXci, sdYci)

        names(sdXYci) <- c("sdX (", "sdX )", "sdY (", "sdY )")
        return(sdXYci)
    }

    ## sd and sd CIs as separate lists for unit of measurement MOA, SMOA, milrad
    sdXY       <- lapply(xyL, function(x) sqrt(diag(cov(x)))) # standard deviations
    sdXYci     <- lapply(xyL, getSDxyCI) # confidence intervals
    res$sdXY   <- Map(makeMOA, sdXY,   dst=dstTarget, conversion=conversion)
    res$sdXYci <- Map(makeMOA, sdXYci, dst=dstTarget, conversion=conversion)

    ## mean distances to group center
    meanDstCtr    <- lapply(xyL, function(x) mean(getDistToCtr(x)))
    meanDstCtrMOA <- Map(makeMOA, meanDstCtr, dst=dstTarget, conversion=conversion)
    res$meanDistToCtr <- do.call("cbind", meanDstCtrMOA)

    ## maximum pairwise distance = maximum group spread
    maxPD      <- lapply(xyL, getMaxPairDist)   # max pairwise distance
    maxSpread  <- lapply(maxPD, function(x) { x$d } )
    maxPDidx   <- sapply(maxPD, function(x) { x$idx } )
    maxSpreadL <- Map(makeMOA, maxSpread, dst=dstTarget, conversion=conversion)
    res$maxPairDist <- do.call("cbind", maxSpreadL)

    ## bounding box figure of merit and diagonal
    ## bbs     <- lapply(xyL, getBoundingBox)   # bounding boxes
    bbs     <- lapply(xyL, getMinBBox)
    bbFoM   <- lapply(bbs, function(x) { x$FoM } )
    bbDiag  <- lapply(bbs, function(x) { x$diag } )
    bbFoML  <- Map(makeMOA, bbFoM,  dst=dstTarget, conversion=conversion)
    bbDiagL <- Map(makeMOA, bbDiag, dst=dstTarget, conversion=conversion)
    res$bbFoM  <- do.call("cbind", bbFoML)
    res$bbDiag <- do.call("cbind", bbDiagL)

    ## for axis limits
    bbMinPts  <- do.call("rbind", lapply(bbs, function(x) x$pts))
    axisLimsX <- c(axisLimsX, bbMinPts[ , 1])
    axisLimsY <- c(axisLimsY, bbMinPts[ , 2])

    ## radius of minimum enclosing circle
    minCircs    <- lapply(xyL, getMinCircle)
    minCircRad  <- lapply(minCircs, function(x) { x$rad } )     # radius
    minCircRadL <- Map(makeMOA, minCircRad, dst=dstTarget, conversion=conversion)
    res$minCircleRad <- do.call("cbind", minCircRadL)

    ## for axis limits
    getMinCircLims <- function(x) {
        cbind(X=c(x$ctr[1] + x$rad, x$ctr[1] - x$rad),
              Y=c(x$ctr[2] + x$rad, x$ctr[2] - x$rad))
    }
    minCircLims <- do.call("rbind", lapply(minCircs, getMinCircLims))
    axisLimsX   <- c(axisLimsX, minCircLims[ , 1])
    axisLimsY   <- c(axisLimsY, minCircLims[ , 2])

    ## 50% circular error probable
    CEPlist <- Map(getCEP, xyL, level=0.5, dstTarget=dstTarget,
                   conversion=conversion, type=CEPtype, accuracy=FALSE)
    CEPl    <- lapply(CEPlist, function(x) { x$CEP[ , CEPtype, drop=FALSE] })
    CEPmat  <- do.call("cbind", CEPl)    # as matrix
    colnames(CEPmat) <- names(xyL)
    res$CEP <- CEPmat

    #####-----------------------------------------------------------------------
    ## tests for equal spread
    ## 2 groups:   Ansari-Bradley for x- and y-coords
    ##             Kruskal-Wallis for distance to center
    ## > 2 groups: Fligner-Killeen for x- and y-coords
    ##             Wilcoxon Rank Sum (= Mann-Whitney U) for distance to center
    dstCtrL   <- lapply(xyL, getDistToCtr) # distances to group center
    dstCtrGrp <- unlist(dstCtrL)           # grouped by Series
    names(dstCtrGrp) <- NULL

    ## create data frame with added Series factor
    dstCtrDf  <- data.frame(dstCtr=dstCtrGrp,
                            Series=factor(rep(1:nlevels(DF$Series), nObs), labels=levels(DF$Series)))

    if(nS == 2) {                        # compare two groups
        res$AnsariX  <- coin::ansari_test(X ~ Series, alternative=ABalt,
                                          data=DF, distribution="exact")
        res$AnsariY  <- coin::ansari_test(Y ~ Series, alternative=ABalt,
                                          data=DF, distribution="exact")
        res$Wilcoxon <- coin::wilcox_test(dstCtr ~ Series, alternative=Walt,
                                          data=dstCtrDf, distribution="exact")
    } else {                             # compare more than two groups
        res$FlignerX <- coin::fligner_test(X ~ Series, data=DF,
                                           distribution=coin::approximate(B=9999))  # x
        res$FlignerY <- coin::fligner_test(Y ~ Series, data=DF,
                                           distribution=coin::approximate(B=9999))  # y
        res$Kruskal  <- coin::kruskal_test(dstCtr ~ Series,    # dist to center
                                           data=dstCtrDf, distribution=coin::approximate(B=9999))
    }

    if(plots) {
        ## infer (x,y)-coord units from conversion
        unitXY  <- getUnits(conversion, first=FALSE)
        unitDst <- getUnits(conversion, first=TRUE)
        devNew  <- getDevice()           # platform-dependent window open

        ## determine axis limits
        xLims <- range(c(DF$X, axisLimsX))
        yLims <- range(c(DF$Y, axisLimsY))

        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        syms <- c(4, 16, 2, 1, 6, 8, 3, 5, 7, 9:13, 15, 17:25)  # data symbols
        cols <- getColors(nS)            # colors

        if(nS > length(syms)) {
            stop(paste("At most", length(syms), "series possible"))
        }

        devNew()                         # open new diagram
        plot(Y ~ X, data=DF, xlim=xLims, ylim=yLims, asp=1, lwd=2,
             pch=syms[unclass(DF$Series)], col=cols[unclass(DF$Series)],
             main="Groups with error ellipses",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim

        ## add confidence ellipses and group centers
        covXY <- lapply(xyL, cov)
        for(i in seq(along=xyL)) {
            drawEllipse(res$ctr[ , i], covXY[[i]], radius=1, fg=cols[i],
                        lwd=2, pch=syms[i], cex=3)
            points(res$ctr[1, i], res$ctr[2, i], pch=syms[i], col=cols[i],
                   cex=3, lwd=2)
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[1:nS],
               col=cols[1:nS], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))

        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(Y ~ X, data=DF, asp=1, xlim=xLims, ylim=yLims, lwd=2,
             pch=syms[unclass(DF$Series)], col=cols[unclass(DF$Series)],
             main="Groups w/ minimum bounding box & maximum spread",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim
        points(res$ctr[1, ], res$ctr[2, ], col=cols[1:nS],
               pch=syms[1:nS], lwd=2, cex=3)

        ## add bounding box and maximum group spread
        for(i in seq(along=xyL)) {
            bb <- bbs[[i]]
            ## drawBox(bb, fg=cols[i])
            drawBox2(bb, fg=cols[i])
            segments(x0=xyL[[i]][maxPDidx[1, i], 1], y0=xyL[[i]][maxPDidx[1, i], 2],
                     x1=xyL[[i]][maxPDidx[2, i], 1], y1=xyL[[i]][maxPDidx[2, i], 2],
                     col=cols[i])
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[1:nS],
               col=cols[1:nS], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))

        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(Y ~ X, data=DF, asp=1, xlim=xLims, ylim=yLims, lwd=2,
             pch=syms[unclass(DF$Series)], col=cols[unclass(DF$Series)],
             main="Groups w/ minimum enclosing circle and mean dist to center",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim
        points(res$ctr[1, ], res$ctr[2, ], col=cols[1:nS], pch=syms[1:nS],
               lwd=2, cex=3)

        ## add circle with mean distance to center and minimum enclosing circle
        for(i in seq(along=xyL)) {
            drawCircle(res$ctr[ , i], radius=meanDstCtr[[i]], fg=cols[i])
            drawCircle(minCircs[[i]], fg=cols[i])
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[1:nS],
               col=cols[1:nS], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
