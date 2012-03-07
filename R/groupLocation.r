groupLocation <-
function(xy, plots=c("0", "1", "2"), conversion="m2cm", dstTarget=25,
         target=c("BDS", "DSB"), unit=c("cm", "in"), caliber=0.9) {
    if(!is.matrix(xy))       { stop("xy must be a matrix") }
    if(!is.numeric(xy))      { stop("xy must be numeric") }
    if(ncol(xy) != 2)        { stop("xy must have two columns") }
    if(!is.numeric(caliber)) { stop("caliber must be numeric") }
    if(caliber <= 0)         { stop("caliber must be positive") }

    plots <- as.character(as.numeric(plots))
    plots <- match.arg(plots)
    plots <- as.numeric(plots)

    #####-----------------------------------------------------------------------
    ## prepare data
    X   <- xy[ , 1]                      # x-coords
    Y   <- xy[ , 2]                      # y-coords
    res <- vector("list", 0)             # empty list to later collect the results

    #####-----------------------------------------------------------------------
    ## location measures
    # library(robustbase)                  # for covMcd()
    res$ctr     <- colMeans(xy)          # center of joint (x,y)-distribution
    res$ctrRob  <- covMcd(xy)$center     # robust estimation of center
    res$distPOA <- c(unit=sqrt(sum(res$ctr^2)),  # distance to point of aim
                     MOA=getMOA(sqrt(sum(res$ctr^2)), dstTarget, conversion))
    res$distPOArob <- c(unit=sqrt(sum(res$ctrRob^2)),
                        MOA=getMOA(sqrt(sum(res$ctrRob^2)), dstTarget, conversion))

    ## Hotelling's T^2 test for equality of (x,y)-center with point of aim (0,0)
    res$Hotelling <- anova(lm(cbind(X, Y) ~ 1), test="Hotelling-Lawley")

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
    nR   <- 1499                         # number of replications
    getM <- function(x, idx) { mean(x[idx]) }  # mean for one replication
    bsX  <- boot(X, statistic=getM, R=nR)      # bootstrap x-coord means
    bsY  <- boot(Y, statistic=getM, R=nR)      # bootstrap y-coord means
    CIbootX <- boot.ci(bsX, conf=1-alpha, type=c("perc", "bca"))  # x-coord
    CIbootY <- boot.ci(bsY, conf=1-alpha, type=c("perc", "bca"))  # y-coord
    res$CIbootX <- rbind(percentile=CIbootX$percent[4:5], BCa=CIbootX$bca[4:5])
    res$CIbootY <- rbind(percentile=CIbootY$percent[4:5], BCa=CIbootY$bca[4:5])
    colnames(res$CIbootX) <- c("(", ")")
    colnames(res$CIbootY) <- c("(", ")")

    if(plots > 0) {
        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        dev.new()                        # open new diagram
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
        points(res$ctr[1],    res$ctr[2],    col="blue", pch=4, lwd=2, cex=1.5)  # group center
        points(res$ctrRob[1], res$ctrRob[2], col="red",  pch=4, lwd=2, cex=1.5)  # robust center
     
        ## add legend
        legend(x="bottomleft", legend=c("group center", "robust group center"),
               col=c("blue", "red"), pch=4, lty=NA, lwd=2, bg="white")
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
