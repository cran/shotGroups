groupLocation <-
function(xy, plots=c("0", "1", "2"), dstTarget=25, conversion="m2cm",
         target=c("BDS", "DSB"), unit=c("cm", "in"), caliber=0.9) {
    if(!is.matrix(xy))       { stop("xy must be a matrix") }
    if(!is.numeric(xy))      { stop("xy must be numeric") }
    if(ncol(xy) != 2)        { stop("xy must have two columns") }
    if(!is.numeric(caliber)) { stop("caliber must be numeric") }
    if(caliber <= 0)         { stop("caliber must be positive") }

    plots <- as.character(as.numeric(plots))
    plots <- match.arg(plots)
    plots <- as.numeric(plots)

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
    ## location measures
    res$ctr <- colMeans(xy)              # center of joint (x,y)-distribution

    ## robust estimation of center
    # library(robustbase)                  # for covMcd()
    if(haveRob) { res$ctrRob <- covMcd(xy)$center }

    res$distPOA <- c(unit=       sqrt(sum(res$ctr^2)),  # distance to point of aim
                      MOA=getMOA(sqrt(sum(res$ctr^2)), dstTarget, conversion))

    if(haveRob) {
        res$distPOArob <- c(unit=       sqrt(sum(res$ctrRob^2)),
                             MOA=getMOA(sqrt(sum(res$ctrRob^2)), dstTarget, conversion))
    }

    ## Hotelling's T^2 test for equality of (x,y)-center with point of aim (0,0)
    if(nrow(xy) > 2) {
        res$Hotelling <- anova(lm(cbind(X, Y) ~ 1), test="Hotelling-Lawley")
    } else { warning("we need >= 3 points for Hotelling T^2") }

    #####-----------------------------------------------------------------------
    ## confidence intervals for x- and y-coords
    ## parametric: t-CI
    alpha <- 0.05                        # alpha-level
    N     <- nrow(xy)                    # number of observations
    tCrit <- qt(c(alpha/2, 1-alpha/2), N-1)  # critical t-values left and right
    Mx    <- mean(X)                     # mean x-coords
    My    <- mean(Y)                     # mean y-coords
    sMx   <- sd(X) / sqrt(N)             # standard error of the mean x
    sMy   <- sd(Y) / sqrt(N)             # standard error of the mean y
    res$CItX <- rev(Mx - tCrit*sMx)      # t-CI x-coords
    res$CItY <- rev(My - tCrit*sMy)      # t-CI y-coords

    ## nonparametric: bootstrap-CIs (percentile and BCa)
    # library(boot)                        # for boot(), boot.ci()
    nR      <- 1499                               # number of replications
    getMean <- function(x, idx) { mean(x[idx]) }  # mean for one replication
    bsX     <- boot(X, statistic=getMean, R=nR)   # bootstrap x-coord means
    bsY     <- boot(Y, statistic=getMean, R=nR)   # bootstrap y-coord means
    CIbootX <- boot.ci(bsX, conf=1-alpha, type=c("perc", "bca"))  # x-coords
    CIbootY <- boot.ci(bsY, conf=1-alpha, type=c("perc", "bca"))  # y-coords
    res$CIbootX <- rbind(percentile=CIbootX$percent[4:5], BCa=CIbootX$bca[4:5])
    res$CIbootY <- rbind(percentile=CIbootY$percent[4:5], BCa=CIbootY$bca[4:5])
    colnames(res$CIbootX) <- c("(", ")")
    colnames(res$CIbootY) <- c("(", ")")

    if(plots > 0) {
        devNew <- getDevice()            # platform-dependent window open
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        if(plots > 1) {                  # show target in background?
            plot(Y ~ X, asp=1, type="n", main="Group (x,y)-coordinates")
            drawTarget(target=target, unit=unit, cex=1.5)
            symbols(Y ~ X, asp=1, main="(x,y)-coordinates", add=TRUE,
                    circles=rep(caliber/2, nrow(xy)), inches=FALSE, 
                    fg=rgb(0.3, 0.3, 0.3, 0.7), bg=rgb(1, 1, 1, 0.5))
        } else {
            plot(Y ~ X, asp=1, main="Group (x,y)-coordinates", pch=16)
            abline(v=0, h=0, col="lightgray")  # add point of aim
        }

        ## add (robust) group center
        points(res$ctr[1], res$ctr[2], col="blue", pch=4, lwd=2, cex=1.5)

        if(haveRob) {
            points(res$ctrRob[1], res$ctrRob[2], col="red",  pch=4, lwd=2, cex=1.5)
        }
     
        ## add legend
        legend(x="bottomleft", legend=c("group center", "robust group center"),
               col=c("blue", "red"), pch=4, lty=NA, lwd=2, bg="white")
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
