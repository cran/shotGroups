groupLocation <-
function(xy, level=0.95, plots=c("0", "1", "2"),
         dstTarget=100, conversion="m2cm",
         target=c("BDS25m", "DSB25m"), caliber=9) {
    UseMethod("groupLocation")
}

groupLocation.data.frame <-
function(xy, level=0.95, plots=c("0", "1", "2"),
         dstTarget=100, conversion="m2cm",
         target=c("BDS25m", "DSB25m"), caliber=9) {
    xy <- getXYmat(xy)
    NextMethod("groupLocation")
}

groupLocation.default <-
function(xy, level=0.95, plots=c("0", "1", "2"),
         dstTarget=100, conversion="m2cm",
         target=c("BDS25m", "DSB25m"), caliber=9) {
    if(!is.matrix(xy))       { stop("xy must be a matrix") }
    if(!is.numeric(xy))      { stop("xy must be numeric") }
    if(ncol(xy) != 2)        { stop("xy must have two columns") }
    if(!is.numeric(level))   { stop("level must be numeric") }
    if(level <= 0)           { stop("level must be > 0") }
    if(!is.numeric(caliber)) { stop("caliber must be numeric") }
    if(caliber <= 0)         { stop("caliber must > 0") }

    plots <- as.character(as.numeric(plots))
    plots <- match.arg(plots)
    plots <- as.numeric(plots)

    ## check if CI level is given in percent
    if(level > 1) {
        while(level > 1) {
            level <- level / 100
        }
        warning(c("level must be in (0,1) and was set to ", level, "\n"))
    }
    
    haveRob <- TRUE                      # can we do robust estimation?
    if(nrow(xy) < 4) {
        warning("We need >= 4 points for robust estimations")
        haveRob <- FALSE
    }                                    # if(nrow(xy) < 4)

    ## infer (x,y)-coord units from conversion and convert caliber size
    unit    <- getUnit(conversion)
    convFac <- getConvFac(paste0("mm2", unit))
    calSize <- convFac * caliber/2

    #####-----------------------------------------------------------------------
    ## prepare data
    X   <- xy[ , 1]                      # x-coords
    Y   <- xy[ , 2]                      # y-coords
    res <- vector("list", 0)             # empty list to later collect the results

    #####-----------------------------------------------------------------------
    ## location measures
    res$ctr <- colMeans(xy)              # center of joint (x,y)-distribution

    ## robust estimation of center
    if(haveRob) {
        res$ctrRob <- robustbase::covMcd(xy)$center
    } else {
        res$ctrRob <- NULL
    }                                    # if(haveRob)

    distPOA <- sqrt(sum(res$ctr^2))      # distance to point of aim
    res$distPOA <- c(unit=distPOA,
                      MOA=getMOA(distPOA, dst=dstTarget, conversion=conversion))

    if(haveRob) {                        # rob distance to point of aim
        distPOArob <- sqrt(sum(res$ctrRob^2))
        res$distPOArob <- c(unit=distPOArob,
                             MOA=getMOA(distPOArob, dst=dstTarget, conversion=conversion))
    } else {
        res$distPOArob <- NULL
    }                                   # if(haveRob)
    
    ## Hotelling's T^2 test for equality of (x,y)-center with point of aim (0,0)
    if(nrow(xy) > 2) {
        res$Hotelling <- anova(lm(cbind(X, Y) ~ 1), test="Hotelling-Lawley")
    } else {
        res$Hotelling <- NULL
        warning("We need >= 3 points for Hotelling's T^2 test")
    }                                    # if(nrow(xy) > 2)

    #####-----------------------------------------------------------------------
    ## confidence intervals for x- and y-coords
    ## parametric: t-CI
    alpha <- 1-level                     # alpha-level
    N     <- nrow(xy)                    # number of observations
    tCrit <- qt(c(alpha/2, 1-alpha/2), N-1)  # critical t-values left and right
    Mx    <- mean(X)                     # mean x-coords
    My    <- mean(Y)                     # mean y-coords
    sMx   <- sd(X) / sqrt(N)             # standard error of the mean x
    sMy   <- sd(Y) / sqrt(N)             # standard error of the mean y
    xCIt  <- rev(Mx - tCrit*sMx)         # t-CI x-coords
    yCIt  <- rev(My - tCrit*sMy)         # t-CI y-coords

    ## nonparametric: bootstrap-CIs for center (basic and BCa)
    nR      <- 1499                                    # number of replications
    getCtr  <- function(x, idx) { colMeans(x[idx, ]) } # center for one replication
    bs      <- boot::boot(xy, statistic=getCtr, R=nR)  # bootstrap centers
    xCIboot <- boot::boot.ci(bs, conf=1-alpha, type=c("basic", "bca"), index=1) # x
    yCIboot <- boot::boot.ci(bs, conf=1-alpha, type=c("basic", "bca"), index=2) # y
    
    ## put both CIs into one component
    res$ctrXci <- rbind(t=xCIt,
                        basic=xCIboot$basic[4:5],
                          BCa=yCIboot$bca[4:5])
    res$ctrYci <- rbind(t=yCIt,
                        basic=yCIboot$basic[4:5],
                          BCa=yCIboot$bca[4:5])
    colnames(res$ctrXci) <- c("x (", "x )")
    colnames(res$ctrYci) <- c("y (", "y )")
    
    if(plots > 0) {
        devNew <- getDevice()            # platform-dependent window open
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        if(plots > 1) {                  # show target in background?
            plot(Y ~ X, asp=1, type="n", main="Group (x,y)-coordinates")
            drawTarget(target, unit=unit, add=TRUE, cex=1.5)
            symbols(Y ~ X, asp=1, main="(x,y)-coordinates", add=TRUE,
                    circles=rep(calSize, nrow(xy)), inches=FALSE, 
                    fg=rgb(0.3, 0.3, 0.3, 0.7), bg=rgb(1, 1, 1, 0.5))
        } else {
            plot(Y ~ X, asp=1, main="Group (x,y)-coordinates", pch=16)
            abline(v=0, h=0, col="lightgray")  # add point of aim
        }                                # if(plots > 1)

        ## add (robust) group center
        points(res$ctr[1], res$ctr[2], col="blue", pch=4, lwd=2, cex=1.5)

        if(haveRob) {
            points(res$ctrRob[1], res$ctrRob[2], col="red",
                   pch=4, lwd=2, cex=1.5)
        }                                # if(haveRob)
     
        ## add legend
        legend(x="bottomleft", legend=c("group center", "robust group center"),
               col=c("blue", "red"), pch=4, lty=NA, lwd=2, bg="white")
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
