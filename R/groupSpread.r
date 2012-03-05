groupSpread <-
function(xy, plots=TRUE, conversion="m2cm", dstTarget=25) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }

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

    ## robust standard deviations of x- and y-coords
    # library(robustbase)                  # for covMcd()
    rob         <- covMcd(xy, cor=TRUE)
    res$sdXYrob <- rbind(unit=sqrt(diag(rob$cov)),   # robust sd
                         MOA=getMOA(sqrt(diag(rob$cov)), dstTarget, conversion))
    colnames(res$sdXYrob) <- c("X", "Y")

    ## covariance- and correlation-matrices, and their robust estimates
    res$covXY    <- cov(xy)              # covariance matrix
    res$covXYrob <- rob$cov              # robust estimate covariance matrix
    res$corXY    <- cor(xy)              # correlation-matrix
    res$corXYrob <- rob$cor              # robust estimate correlation-matrix

    ## more spread measures
    dstCtr  <- getDistToCtr(xy)          # distances to group center
    maxPD   <- getMaxPairDist(xy)        # maximum group spread
    # bb      <- getBoundingBox(xy)        # bounding box
    bb      <- getMinBBox(xy)            # minimum bounding box
    minCirc <- getMinCircle(xy)          # minimum enclosing circle
    res$meanDistToCtr <- c(unit=mean(dstCtr),  # mean distance to center
                           MOA=getMOA(mean(dstCtr), dstTarget, conversion))
    res$maxPairDist   <- c(unit=maxPD$d,       # maximum pairwise distance
                           MOA=getMOA(maxPD$d, dstTarget, conversion))
    res$groupWidth    <- c(unit=bb$width,      # group width
                           MOA=getMOA(bb$width, dstTarget, conversion))
    res$groupHeight   <- c(unit=bb$height,     # group height
                           MOA=getMOA(bb$height, dstTarget, conversion))
    res$minCircleRad  <- c(unit=minCirc$rad,   # radius min enclosing circle
                           MOA=getMOA(minCirc$rad, dstTarget, conversion))

    if(plots) {
        #####-------------------------------------------------------------------
        ## diagram: histogram for distances to group center
        dev.new()                        # open new diagram
        hist(dstCtr, breaks="FD", freq=FALSE,
             main="Histogram distances to center w/ kernel density estimate")
        rug(jitter(dstCtr))              # show single values
        lines(density(dstCtr), lwd=2, col="red")  # add kernel density estimate
        legend(x="topleft", legend="kernel density estimate", col="red",
               lty=1, lwd=2, bg="white")
        
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        ## determine axis limits
        xLims <- range(c(bb$pts[ , 1], minCirc$ctr[1] + c(-minCirc$rad, minCirc$rad)))
        yLims <- range(c(bb$pts[ , 2], minCirc$ctr[2] + c(-minCirc$rad, minCirc$rad)))

        ## calculate radius for confidence ellipse
        alpha <- 0.05                    # alpha-level
        N     <- nrow(xy)                # number of observations
        dfn   <- ncol(xy)                # numerator df
        dfd   <- N-1                     # denominator df
        rad   <- sqrt(dfn*qf(1-alpha, dfn, dfd))  # radius
        dev.new()                        # open new diagram
        plot(Y ~ X, asp=1, xlim=xLims, ylim=yLims, pch=20, main="Group (x,y)-coordinates")
        abline(v=0, h=0, col="lightgray")  # add point of aim

        ## add group center and robust estimate for group center
        ctr      <- colMeans(xy)         # group center
        covXY    <- cov(xy)              # covariance matrix (x,y)-coords
        ctrRob   <- rob$center           # robust estimate: group center
        covXYrob <- rob$cov              # robust estimate: group covariance matrix

        points(ctr[1],    ctr[2],    col="red",  pch=4, lwd=2, cex=1.5)
        points(ctrRob[1], ctrRob[2], col="blue", pch=4, lwd=2, cex=1.5)
        
        ## add confidence ellipses (parametric, robust),
        ## and a circle with mean distance to center
        drawEllipse(ctr,    covXY,    radius=rad, pch=4, fg="red",  lwd=2)
        drawEllipse(ctrRob, covXYrob, radius=rad, pch=4, fg="blue", lwd=2)
        drawCircle(ctr[1], ctr[2], radius=mean(dstCtr), fg="black", lwd=2)

        ## add legend
        legend(x="bottomleft", legend=c("center", "center (robust)",
               "95% confidence ellipse", "95% confidence ellipse (robust)",
               "mean distance to center"),
               col=c("red", "blue", "red", "blue", "black"), pch=c(4, 4, NA, NA, NA),
               lty=c(NA, NA, 1, 1, 1), lwd=2, bg="white")
        
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        dev.new()                        # open new diagram
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
