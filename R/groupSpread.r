groupSpread <-
function(xy, plots=TRUE, level=0.5, dstTarget=25, conversion="m2cm") {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }

    haveRob <- TRUE                      # can we do robust estimation?
    if(nrow(xy) < 4) {
        warning("we need >= 4 points for robust estimations")
        haveRob <- FALSE
    }

    #####-----------------------------------------------------------------------
    ## prepare data
    X   <- xy[ , 1]                      # x-coords
    Y   <- xy[ , 2]                      # y-coords
    res <- vector("list", 0)             # empty list to later collect the results

    #####-----------------------------------------------------------------------
    ## regular standard deviations of x- and y-coords
    res$sdXY <- rbind(unit=c(sd(X), sd(Y)),      # standard deviations
                      MOA=getMOA(c(sd(X), sd(Y)), dstTarget, conversion))
    colnames(res$sdXY) <- c("X", "Y")

    ## and their 95%-confidence intervals
    N     <- nrow(xy)
    alpha <- 0.05
    xLow  <- sqrt((N-1)*var(X) / qchisq(1-(alpha/2), N-1))
    xHi   <- sqrt((N-1)*var(X) / qchisq(   alpha/2,  N-1))
    yLow  <- sqrt((N-1)*var(Y) / qchisq(1-(alpha/2), N-1))
    yHi   <- sqrt((N-1)*var(Y) / qchisq(   alpha/2,  N-1))
    res$CIsdXY <- rbind(sdX=c("("=xLow, ")"=xHi), sdY=c("("=yLow, ")"=yHi))

    ## robust standard deviations of x- and y-coords
    if(haveRob) {
        # library(robustbase)              # for covMcd()
        rob         <- covMcd(xy, cor=FALSE)
        res$sdXYrob <- rbind(unit=sqrt(diag(rob$cov)),   # robust sd
                              MOA=getMOA(sqrt(diag(rob$cov)), dstTarget, conversion))
        colnames(res$sdXYrob) <- c("X", "Y")
    }

    ## (robust) center and covariance-matrix
    ctr <- colMeans(xy)                  # group center
    if(haveRob) { ctrRob <- rob$center }
    res$covXY <- cov(xy)                 # covariance matrix
    if(haveRob) { res$covXYrob <- rob$cov }

    ## mean distance to group center
    dstCtr <- getDistToCtr(xy)
    res$meanDistToCtr <- c(unit=mean(dstCtr),
                           MOA=getMOA(mean(dstCtr), dstTarget, conversion))

    ## maximum pairwise distance
    maxPD <- getMaxPairDist(xy)
    res$maxPairDist <- c(unit=maxPD$d, MOA=getMOA(maxPD$d, dstTarget, conversion))

    ## width and height of minimum bounding box
    # bb            <- getBoundingBox(xy)   # bounding box
    bb            <- getMinBBox(xy)      # minimum bounding box
    groupWidth    <- c(bb$width,  getMOA(bb$width,  dstTarget, conversion))
    groupHeight   <- c(bb$height, getMOA(bb$height, dstTarget, conversion))
    res$groupRect <- cbind(width=groupWidth, height=groupHeight)
    rownames(res$groupRect) <- c("unit", "MOA")

    ## radius of minimum enclosing circle
    minCirc <- getMinCircle(xy)          # minimum enclosing circle
    res$minCircleRad <- c(unit=minCirc$rad,
                          MOA=getMOA(minCirc$rad, dstTarget, conversion))

    ## confidence ellipse measures
    confEll <- getConfEll(xy, level, dstTarget, conversion)
    res$confEll <- confEll$size
    if(haveRob) { res$confEllRob <- confEll$sizeRob }
    res$confEllShape <- confEll$shape
    if(haveRob) { res$confEllShapeRob <- confEll$shapeRob }

    ## circular error probable
    res$CEPrand <- getCEP(xy, dstTarget, conversion)$RAND

    if(plots) {
        devNew <- getDevice()            # platform-dependent window open
        #####-------------------------------------------------------------------
        ## diagram: histogram for distances to group center
        devNew()                         # open new diagram
        hist(dstCtr, breaks="FD", freq=FALSE,
             main="Histogram distances to center w/ kernel density estimate")
        rug(jitter(dstCtr))              # show single values
        lines(density(dstCtr), lwd=2, col="red")  # add kernel density estimate
        legend(x="topright", legend="kernel density estimate", col="red",
               lty=1, lwd=2, bg="white")
        
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        ## determine axis limits
        xLims <- range(c(bb$pts[ , 1], minCirc$ctr[1] + c(-minCirc$rad, minCirc$rad)))
        yLims <- range(c(bb$pts[ , 2], minCirc$ctr[2] + c(-minCirc$rad, minCirc$rad)))

        devNew()                         # open new diagram
        plot(Y ~ X, asp=1, xlim=xLims, ylim=yLims, pch=20, main="Group (x,y)-coordinates")
        abline(v=0, h=0, col="lightgray")  # add point of aim

        ## add group center and robust estimate for group center
        points(ctr[1], ctr[2], col="red",  pch=4, lwd=2, cex=1.5)
        if(haveRob) { points(ctrRob[1], ctrRob[2], col="blue", pch=4, lwd=2, cex=1.5) }
        
        ## add confidence ellipses (parametric, robust),
        ## and a circle with mean distance to center
        drawEllipse(ctr, res$covXY, radius=confEll$magFac, pch=4, fg="red",  lwd=2)
        if(haveRob) {
            drawEllipse(ctrRob, res$covXYrob, radius=confEll$magFac, pch=4, fg="blue", lwd=2)
        }
        drawCircle(ctr[1], ctr[2], radius=mean(dstCtr), fg="black", lwd=2)

        ## add legend
        legend(x="bottomleft", legend=c("center", "center (robust)",
               paste(100*level, "% confidence ellipse", sep=""),
               paste(100*level, "% confidence ellipse (robust)", sep=""),
               "mean distance to center"),
               col=c("red", "blue", "red", "blue", "black"), pch=c(4, 4, NA, NA, NA),
               lty=c(NA, NA, 1, 1, 1), lwd=2, bg="white")
        
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(Y ~ X, asp=1, xlim=xLims, ylim=yLims, pch=20, main="Group (x,y)-coordinates")
        abline(v=0, h=0, col="lightgray")                             # add point of aim
        points(ctr[1], ctr[2], col="magenta", pch=4, lwd=2, cex=1.5)  # add group center

        ## add bounding box, minimum enclosing circle, and maximum group spread
        ## drawBox(bb$pts[1], bb$pts[2], bb$pts[3], bb$pts[4], fg="red", lwd=2)
        drawBox2(bb$pts, fg="red", lwd=2)
        drawCircle(minCirc$ctr[1], minCirc$ctr[2], minCirc$rad, fg="blue", lwd=2)
        segments(x0=xy[maxPD$idx[1], 1], y0=xy[maxPD$idx[1], 2],
                 x1=xy[maxPD$idx[2], 1], y1=xy[maxPD$idx[2], 2], col="green3", lwd=2)
        
        ## add legend
        legend(x="bottomleft", legend=c("group center", "minimum bounding box",
               "minimum enclosing circle", "maximum group spread"),
               col=c("magenta", "red", "blue", "green3"),
               pch=c(4, NA, NA, NA), lty=c(NA, 1, 1, 1), lwd=2, bg="white")
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
