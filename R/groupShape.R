groupShape <-
function(xy, plots=TRUE, bandW=0.5, outlier=c("mcd", "pca"), ...) {
    UseMethod("groupShape")
}

groupShape.data.frame <-
function(xy, plots=TRUE, bandW=0.5, outlier=c("mcd", "pca"), ...) {
    xy <- getXYmat(xy)
    NextMethod("groupShape")
}

groupShape.default <-
function(xy, plots=TRUE, bandW=0.5, outlier=c("mcd", "pca"), ...) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }

    outlier <- match.arg(outlier)

    haveRob <- TRUE                      # can we do robust estimation?
    if(nrow(xy) < 4) {
        warning(c("We need >= 4 points for robust estimations,\n",
                  "outlier analysis, and chi^2 plot for multivariate normality"))
        haveRob <- FALSE
    }                                    # if(nrow(xy) < 4)
    
    #####-----------------------------------------------------------------------
    ## prepare data
    X   <- xy[ , 1]                      # x-coords
    Y   <- xy[ , 2]                      # y-coords
    res <- vector("list", 0)             # empty list to later collect the results

    if(plots) {
        devNew <- getDevice()            # platform-dependent window open
    }                                    # if(plots)

    ## correlation matrix of (x,y)-coords
    res$corXY <- cor(xy)

    if(haveRob) {
        rob <- robustbase::covMcd(xy, cor=TRUE)
        res$corXYrob <- rob$cor          # robust estimate correlation-matrix

        #####-------------------------------------------------------------------
        ## outlier-analysis for joint distribution of (x,y)-coords
        devNew()                         # open new diagram
        op <- par(no.readonly=TRUE)      # save device parameters
        outXY <- mvoutlier::aq.plot(xy)  # outlier-analysis-plot
        res$Outliers <- if(outlier == "mcd") {  #  identified outliers
            which(outXY$outliers)
        } else if(outlier == "pca") {
            which(mvoutlier::pcout(xy, makeplot=FALSE, ...)$wfinal01 == 0)
        }
        par(op)                          # reset device parameters
    } else {
        res$corXYrob <- NULL
        res$Outliers <- NULL
    }                                    # if(haveRob)
    
    #####-----------------------------------------------------------------------
    ## normality tests
    ## Shapiro-Wilk-Tests for normality of (x,y)-coords separately
    if(nrow(xy) >= 3) {
        res$ShapiroX <- shapiro.test(X)  # normality x-coords
        res$ShapiroY <- shapiro.test(Y)  # normality y-coords

        ## Energy-Test for multivariate normality of joint (x,y)-distribution
        res$multNorm <- energy::mvnorm.etest(xy)
    } else {
        res$ShapiroX <- NULL
        res$ShapiroY <- NULL
        res$multNorm <- NULL
        warning("We need >= 3 points for normality tests")
    }                                    # if(nrow(xy) >= 3)

    if(plots) {
        #####-------------------------------------------------------------------
        ## diagram: separate Q-Q-plots for eyeballing normality in x- and y-coords
        devNew()                         # open new diagram
        qqnorm(X, pch=20, main="Q-Q-plot x-coordinates for eyeballing normality",
               xlab="Quantiles from normal distribution")
        qqline(X, col="red", lwd=2)      # reference line
        devNew()                         # open new diagram
        qqnorm(Y, pch=20, main="Q-Q-plot y-coordinates for eyeballing normality",
               xlab="Quantiles from normal distribution")
        qqline(Y, col="red", lwd=2)      # reference line
        
        #####-------------------------------------------------------------------
        ## diagram: histograms for x- and y-coords
        ## x-coords
        ## choose y-axis limits
        maxNorm <- getMaxNorm(X, 2)[2]
        dens    <- hist(X, breaks="FD", plot=FALSE)$density
        yLims   <- c(0, max(c(dens, maxNorm)))

        devNew()                         # open new diagram
        hist(X, ylim=yLims, breaks="FD", freq=FALSE,
             main="Histogram x-coordinates w/ kernel density estimate")
        rug(jitter(X))                   # show single values
        
        ## add fitted normal curve and kernel density estimate
        dnormX <- function(x) { dnorm(x, mean(X), sd(X)) }
        curve(dnormX, lwd=2, col="blue", add=TRUE)
        lines(density(X), lwd=2, col="red")  # kernel density estimate
        
        ## add legend
        legend(x="topleft", legend=c("normal distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=c(1, 1), lwd=c(2, 2), bg="white")
        
        ## histogram y-coords
        ## choose y-axis limits
        maxNorm <- getMaxNorm(Y, 2)[2]
        dens    <- hist(Y, breaks="FD", plot=FALSE)$density
        yLims   <- c(0, max(c(dens, maxNorm)))

        devNew()                         # open new diagram
        hist(Y, ylim=yLims, breaks="FD", freq=FALSE,
             main="Histogram y-coordinates w/ kernel density estimate")
        rug(jitter(Y))                   # show single values
        
        ## add fitted normal curve and kernel density estimate
        dnormY <- function(x) { dnorm(x, mean(Y), sd(Y)) }
        curve(dnormY, lwd=2, col="blue", add=TRUE)
        lines(density(Y), lwd=2, col="red")  # kernel density estimate
        
        ## add legend
        legend(x="topleft", 
               legend=c("normal distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=c(1, 1), lwd=c(2, 2), bg="white")

        ## chi-square qq-plot for eyeballing multivariate normality
        ## quantiles of squared robust Mahalanobis distance against quantiles
        ## from chi^2 distribution with 2 df
        if(haveRob) {
            ctrRob   <- rob$center       # robust estimate: group center,
            covXYrob <- rob$cov          # group covariance matrix

            ## squared robust Mahalanobis-distance
            mDstSq <- mahalanobis(xy, center=ctrRob, cov=covXYrob)
            plot(qchisq(ppoints(mDstSq), df=ncol(xy)), sort(mDstSq),
                 xlab="Quantiles chi^2 distribution",
                 ylab="Quantiles (robust Mahalanobis distances)^2", pch=20,
                 main="Chi^2 Q-Q-plot for eyeballing multivariate normality")
            abline(a=0, b=1, col="red", lwd=2)  # add a reference line
        }                                # if(haveRob)

        #####-------------------------------------------------------------------
        ## diagram: 2D-kernel density estimate for joint (x,y)-distribution
        devNew()                         # open new diagram
        smoothScatter(X, Y, asp=1, bandwidth=bandW,
                      main="2D-kernel density estimate and error ellipses")
        abline(h=0, v=0, lwd=2)          # add point of aim

        ## add group center and robust estimate for group center
        if(haveRob) {
            points(ctrRob[1], ctrRob[2], col=rgb(0.5, 0.5, 0.5, 0.4),
                   pch=4, lwd=2, cex=1.5)
        
            ## add robust error ellipses with radius = 1 and = 2
            drawEllipse(ctrRob, covXYrob, radius=1, fg=rgb(0.5, 0.5, 0.5, 0.4),
                        pch=4, lwd=2)
            drawEllipse(ctrRob, covXYrob, radius=2, fg=rgb(0.5, 0.5, 0.5, 0.4),
                        pch=4, lwd=2)
        }                                # if(haveRob)
        
        ## add legend
        legend(x="bottomleft", legend=c("robust center", "robust error ellipses"),
               col=c("darkgray", "darkgray"), pch=c(4, NA), lty=c(NA, 1),
               lwd=c(2, 2), bg="white")
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
