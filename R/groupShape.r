groupShape <-
function(xy, plots=TRUE, bandW=0.5) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }

    #####-----------------------------------------------------------------------
    ## prepare data
    X   <- xy[ , 1]                      # x-coords
    Y   <- xy[ , 2]                      # y-coords
    res <- vector("list", 0)             # empty list to later collect the results

    #####-----------------------------------------------------------------------
    ## outlier-analysis for joint distribution of (x,y)-coords
    dev.new()                            # open new diagram
    # library(mvoutlier)                   # for aq.plot()
    outXY        <- aq.plot(xy)          # outlier-analysis-plot
    res$Outliers <- which(outXY$outliers)  #  identified outliers

    #####-----------------------------------------------------------------------
    ## normality tests
    ## Shapiro-Wilk-Tests for normality of (x,y)-coords separately
    res$ShapiroX <- shapiro.test(X)      # normality x-coords
    res$ShapiroY <- shapiro.test(Y)      # normality y-coords

    ## Energy-Test for multivariate normality of joint (x,y)-distribution
    # library(energy)                    # for mvnorm.etest()
    res$multNorm <- mvnorm.etest(xy)

    if(plots) {
        #####-------------------------------------------------------------------
        ## diagram: separate Q-Q-plots for eyeballing normality in x- and y-coords
        dev.new()                        # open new diagram
        qqnorm(X, pch=20, main="Q-Q-plot x-coordinates for eyeballing normality")  # x-coords
        qqline(X, col="red", lwd=2)      # reference line
        dev.new()                        # open new diagram
        qqnorm(Y, pch=20, main="Q-Q-plot y-coordinates for eyeballing normality")  # y-coords
        qqline(Y, col="red", lwd=2)      # reference line
        
        #####-------------------------------------------------------------------
        ## diagram: histograms for x- and y-coords
        ## x-coords
        yLims <- c(0, getMaxNorm(X, 1)[2])  # choose y-axis limits
        
        dev.new()                        # open new diagram
        hist(X, ylim=yLims, breaks="FD", freq=FALSE,
             main="Histogram x-coordinates w/ kernel density estimate")
        rug(jitter(X))                    # show single values
        
        ## add fitted normal curve and kernel density estimate
        curve(dnorm(x, mean(X), sd(X)), lwd=2, col="blue", add=TRUE)
        lines(density(X), lwd=2, col="red")  # kernel density estimate
        
        ## add legend
        legend(x="topleft", legend=c("normal distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=c(1, 1), lwd=c(2, 2), bg="white")
        
        ## histogram y-coords
        yLims <- c(0, getMaxNorm(Y, 1)[2])  # choose y-axis limits
        dev.new()                        # open new diagram
        hist(Y, ylim=yLims, breaks="FD", freq=FALSE,
             main="Histogram y-coordinates w/ kernel density estimate")
        rug(jitter(Y))                   # show single values
        
        ## add fitted normal curve and kernel density estimate
        curve(dnorm(x, mean(Y), sd(Y)), lwd=2, col="blue", add=TRUE)
        lines(density(Y), lwd=2, col="red")  # kernel density estimate
        
        ## add legend
        legend(x="topleft", legend=c("normal distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=c(1, 1), lwd=c(2, 2), bg="white")

        ## chi-square plot for eyeballing multivariate normality
        ## quantiles of robust squared Mahalanobis distance against quantiles
        ## from chi^2 distribution with 2 df
        # library(robustbase)              # for covMcd()
        rob      <- covMcd(xy)
        ctrRob   <- rob$center           # robust estimate: group center
        covXYrob <- rob$cov              # robust estimate: group covariance matrix

        ## squared robust Mahalanobis-distance
        mDist <- mahalanobis(xy, center=ctrRob, cov=covXYrob)
        sMd   <- sort(mDist)
        qq    <- (0.5:length(mDist)) / length(mDist)
        qChi  <- qchisq(qq, df=ncol(xy))
        dev.new()                        # open new diagram
        plot(sMd, qChi, xlab="Quantiles squared robust Mahalanobis distance",
             ylab="Quantiles of corresponding Chi^2 distrib", 
             main="Chi^2 plot for eyeballing multivariate normality", col=3)

        ## add a reference line throught first and third quantile
        y     <- quantile(mDist, c(0.25, 0.75))
        x     <- qchisq(c(0.25, 0.75), df=ncol(xy))
        slope <- diff(y) / diff(x)
        int   <- y[1] - slope*x[1]
        abline(int, slope, col="black", lwd=2)
        
        #####-------------------------------------------------------------------
        ## diagram: 2D-kernel density estimate for joint (x,y)-distribution
        dev.new()                        # open new diagram
        smoothScatter(X, Y, asp=1, bandwidth=bandW,
                      main="2D-kernel density estimate and characteristic ellipses")
        abline(h=0, v=0, lwd=2)          # add point of aim

        ## add group center and robust estimate for group center
        points(ctrRob[1], ctrRob[2], col=rgb(0.5, 0.5, 0.5, 0.4), pch=4, lwd=2, cex=1.5)
        
        ## add robust characteristic ellipsoid with radius = 1 and = 2
        drawEllipse(ctrRob, covXYrob, radius=1, fg=rgb(0.5, 0.5, 0.5, 0.4), pch=4, lwd=2)
        drawEllipse(ctrRob, covXYrob, radius=2, fg=rgb(0.5, 0.5, 0.5, 0.4), pch=4, lwd=2)
        
        ## add legend
        legend(x="bottomleft", legend=c("robust center", "robust characteristic ellipses"),
               col=c("darkgray", "darkgray"), pch=c(4, NA), lty=c(NA, 1),
               lwd=c(2, 2), bg="white")
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
