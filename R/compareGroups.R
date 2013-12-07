compareGroups <-
function(DF, plots=TRUE, xyTopLeft=TRUE, ABalt=c("two.sided", "less", "greater"),
         Walt=c("two.sided", "less", "greater"), conversion="m2cm") {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    #####-----------------------------------------------------------------------
    ## make sure DF has the required variable names and at least two groups
    varNames <- names(DF)                # what variables are present
    needs    <- c("Series", "Distance", "Aim.X", "Aim.Y", "Point.X", "Point.Y")  # required
    has      <- needs %in% varNames      # what we have of the required ones
    if(!all(has)) {
        stop(cat("the data frame is missing variable(s)\n", needs[!has], "\n"))
    }

    #####-----------------------------------------------------------------------
    ## prepare data
    res <- vector("list", 0)             # empty list to later collect the results
    DF$Series <- droplevels(DF$Series)   # remove all non-used factor levels

    ## check if we have enough groups and points per group
    if(nlevels(DF$Series) < 2)    { stop("we need >= 2 groups for a comparison") }
    if(any(table(DF$Series) < 2)) { stop("we need >= 2 points in each group") }

    ## prepare data: coords relative to point of aim
    ## y-coords exported from OnTarget: (0,0) is top-left
    DF$X <- DF$Point.X - DF$Aim.X        # x-coords
    if(xyTopLeft) {
        DF$Y <- -(DF$Point.Y - DF$Aim.Y)
    } else {
        DF$Y <-   DF$Point.Y - DF$Aim.Y
    }
    xy        <- cbind(DF$X, DF$Y)       # new (x,y)-coords as matrix
    dstTarget <- tapply(DF$Distance, DF$Series, mean)  # distances to target

    ## for each group extract the new (x,y)-coords as a matrix
    xyL <- lapply(split(DF, DF$Series), function(x) {
                  data.matrix(subset(x, select=c("X", "Y"))) } )
    nS <- length(xyL)                    # total number of series
    names(xyL) <- paste("Series", levels(DF$Series), sep="")

    #####-----------------------------------------------------------------------
    ## location measures
    res$ctr     <- sapply(xyL, colMeans)  # group centers
    distPOA     <- sqrt(colSums(res$ctr^2))  # distances to point of aim
    distPOAmoa  <- getMOA(distPOA, dstTarget, conversion)
    res$distPOA <- rbind(unit=distPOA, MOA=distPOAmoa)

    ## multivariate location test for equal group centers
    res$MANOVA <- anova(lm(cbind(X, Y) ~ Series, data=DF), test="Wilks")

    #####-----------------------------------------------------------------------
    ## shape measures
    ## correlation matrices for x- and y-coords
    res$corXY <- lapply(xyL, cor)

    #####-----------------------------------------------------------------------
    ## spread measures
    ## standard deviations for x- and y-coords
    covXY    <- lapply(xyL, cov)         # covariance matrices (x,y)-coords
    res$sdXY <- sapply(covXY, function(x) sqrt(diag(x)))  # std devs

    ## mean distances to group center
    dstList   <- lapply(xyL, getDistToCtr) # distances to group center
    dstCtr    <- do.call("c", dstList)     # as vector
    dstCtrMOA <- getMOA(dstCtr, DF$Distance, conversion)     # as MOA
    meanDistToCtr     <- tapply(dstCtr,    DF$Series, mean)  # mean distance
    meanDistToCtrMOA  <- tapply(dstCtrMOA, DF$Series, mean)  # as MOA
    res$meanDistToCtr <- rbind(unit=meanDistToCtr, MOA=meanDistToCtrMOA)
    colnames(res$meanDistToCtr) <- names(xyL)

    ## maximum pairwise distance = maximum group spread
    maxPD        <- lapply(xyL, getMaxPairDist)  # max pairwise distance
    maxSpread    <- sapply(maxPD, function(x) { x$d } )
    maxPDidx     <- sapply(maxPD, function(x) { x$idx } )
    maxSpreadMOA <- getMOA(maxSpread, dstTarget, conversion)  # as MOA
    res$maxPairDist <- rbind(unit=maxSpread, MOA=maxSpreadMOA)

    ## bounding box width and height
    ## bbs             <- lapply(xyL, getBoundingBox)   # bounding boxes
    bbs             <- lapply(xyL, getMinBBox)
    groupWidth      <- sapply(bbs, function(x) { x$width } )
    groupHeight     <- sapply(bbs, function(x) { x$height } )
    groupWidthMOA   <- getMOA(groupWidth,  dstTarget, conversion)
    groupHeightMOA  <- getMOA(groupHeight, dstTarget, conversion)
    res$groupWidth  <- rbind(unit=groupWidth,  MOA=groupWidthMOA)
    res$groupHeight <- rbind(unit=groupHeight, MOA=groupHeightMOA)

    ## radius of minimum enclosing circle
    minCircs         <- lapply(xyL, getMinCircle)
    minCircleRad     <- sapply(minCircs, function(x) { x$rad } )     # radius
    minCircleRadMOA  <- getMOA(minCircleRad, dstTarget, conversion)  # as MOA
    res$minCircleRad <- rbind(unit=minCircleRad, MOA=minCircleRadMOA)

    ## 50% circular error probable
    CEPlist <- vector("list", nS)
    for(i in seq(along=numeric(nS))) {
        CEPlist[[i]] <- getCEP(xyL[[i]], dstTarget[i], conversion)
    }
    CEPrandL <- lapply(CEPlist, function(x) { x$RAND[ , "50%", drop=FALSE] } )
    CEPmat   <- do.call("cbind", CEPrandL)    # as matrix
    colnames(CEPmat) <- names(xyL)
    res$CEPrand <- CEPmat

    #####-----------------------------------------------------------------------
    ## tests for equal spread
    ## 2 groups:   Ansari-Bradley for x- and y-coords
    ##             Kruskal-Wallis for distance to center
    ## > 2 groups: Fligner-Killeen for x- and y-coords
    ##             Wilcoxon Rank Sum (= Mann-Whitney U) for distance to center
    if(nS == 2) {                        # compare two groups
        res$AnsariX  <- coin::ansari_test(     X ~ Series, alternative=ABalt,
                                          data=DF, distribution="exact")
        res$AnsariY  <- coin::ansari_test(     Y ~ Series, alternative=ABalt,
                                          data=DF, distribution="exact")
        res$Wilcoxon <- coin::wilcox_test(dstCtr ~ Series, alternative=Walt,
                                          data=DF, distribution="exact")
    } else {                             # compare more than two groups
        res$FlignerX <- coin::fligner_test(X ~ Series, data=DF,
                                           distribution="approximate")  # x
        res$FlignerY <- coin::fligner_test(Y ~ Series, data=DF,
                                           distribution="approximate")  # y
        res$Kruskal  <- coin::kruskal_test(dstCtr ~ Series,    # dist to center
                                           data=DF, distribution="approximate")
    }

    if(plots) {
        devNew <- getDevice()            # platform-dependent window open
        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        syms <- c(4, 16, 2, 1, 6, 8, 3, 5, 7, 9:13, 15, 17:25)  # data symbols
        cols <- getColors(nS)            # colors
        
        if(nS > length(syms)) { stop(paste("at most", length(syms), "series possible")) }
        
        devNew()                         # open new diagram
        plot(DF$X, DF$Y, asp=1, lwd=2, xlab="X", ylab="Y",
             pch=syms[unclass(DF$Series)], col=cols[unclass(DF$Series)],
             main="Groups with error ellipses")
        abline(v=0, h=0, col="lightgray")  # add point of aim
        
        ## add confidence ellipses and group centers
        for(i in seq(along=numeric(nS))) {
            drawEllipse(res$ctr[ , i], covXY[[i]], radius=1, fg=cols[i],
                        lwd=2, pch=syms[i], cex=3)
            points(res$ctr[1, i], res$ctr[2, i], pch=syms[i], col=cols[i], cex=3, lwd=2)
        }
        
        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[1:nS],
               col=cols[1:nS], lwd=2, pt.cex=1.5, bg="white")
        
        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        ## determine axis limits
        xLims <- range(c(DF$X, bbs[[1]]$pts[ , 1],
                       minCircs[[1]]$ctr[1] + c(-minCircs[[1]]$rad, minCircs[[1]]$rad)))
        yLims <- range(c(DF$Y, bbs[[1]]$pts[ , 2],
                       minCircs[[1]]$ctr[2] + c(-minCircs[[1]]$rad, minCircs[[1]]$rad)))

        devNew()                         # open new diagram
        plot(DF$X, DF$Y, asp=1, xlim=xLims, ylim=yLims, xlab="X", ylab="Y", lwd=2,
             pch=syms[unclass(DF$Series)], col=cols[unclass(DF$Series)],
             main="Groups w/ minimum bounding box & maximum spread")
        abline(v=0, h=0, col="lightgray")  # add point of aim
        points(res$ctr[1, ], res$ctr[2, ], col=cols[1:nS], pch=syms[1:nS], lwd=2, cex=3)
        
        ## add bounding box and maximum group spread
        for(i in seq(along=numeric(nS))) {
            bb <- bbs[[i]]
            ## drawBox(bb$pts[1], bb$pts[2], bb$pts[3], bb$pts[4], fg=cols[i])
            drawBox2(bb$pts, fg=cols[i])
            segments(x0=xyL[[i]][maxPDidx[1, i], 1], y0=xyL[[i]][maxPDidx[1, i], 2],
                     x1=xyL[[i]][maxPDidx[2, i], 1], y1=xyL[[i]][maxPDidx[2, i], 2],
                     col=cols[i])
        }
        
        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[1:nS],
               col=cols[1:nS], lwd=2, pt.cex=1.5, bg="white")
        
        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(DF$X, DF$Y, asp=1, xlim=xLims, ylim=yLims, xlab="X", ylab="Y", lwd=2,
             pch=syms[unclass(DF$Series)], col=cols[unclass(DF$Series)],
             main="Groups w/ minimum enclosing circle and mean dist to center")
        abline(v=0, h=0, col="lightgray")  # add point of aim
        points(res$ctr[1, ], res$ctr[2, ], col=cols[1:nS], pch=syms[1:nS], lwd=2, cex=3)
        
        ## add circle with mean distance to center and minimum enclosing circle
        for(i in seq(along=numeric(nS))) {
            drawCircle(res$ctr[1, i], res$ctr[2, i], radius=meanDistToCtr[i], fg=cols[i])
            mc <- minCircs[[i]]
            drawCircle(mc$ctr[1], mc$ctr[2], radius=mc$rad, fg=cols[i])
        }
        
        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[1:nS],
               col=cols[1:nS], lwd=2, pt.cex=1.5, bg="white")
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
