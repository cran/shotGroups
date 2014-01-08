getCEP <-
function(xy, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type=c("Rayleigh", "Grubbs", "RAND")) {
    UseMethod("getCEP")
}

getCEP.data.frame <-
function(xy, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type=c("Rayleigh", "Grubbs", "RAND")) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getCEP")
}

getCEP.default <-
function(xy, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type=c("Rayleigh", "Grubbs", "RAND")) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }

    if(!all(type %in% c("Rayleigh", "Grubbs", "RAND"))) {
        stop('type must be any of "Rayleigh", "Grubbs", "RAND"')
    }                                    # if(...)

    ## some basic calculations used later
    level <- 1-0.5                       # 50% CEP
    N     <- nrow(xy)                    # number of points
    ctr   <- colMeans(xy)                # center
    eig   <- eigen(cov(xy))              # eigenvalues and eigenvectors

    ## make sure eigenvalues >= 0 when very small
    eigVal <- eig$values * sign(eig$values)

    ## error ellipse characteristics
    ellRad <- sqrt(eigVal)               # ellipse major radii
    aspRat <- ellRad[1] / ellRad[2]      # aspect ratio of ellipse
    flat   <- 1 - (ellRad[2]/ellRad[1])  # flattening of ellipse

    ## modified RAND-234 CEP estimate for 50% from Williams, 1997
    ## using the semi-major and semi-minor axes of the error ellipse (PCA)
    if(aspRat > 4) {
        warning(c("\nAspect ratio of error ellipse is ",
                  round(aspRat, 2) , " (> 4),\n",
                  "more than what RAND CEP was considered for\n"))
    }

    RAND50b <- 0.563*ellRad[1] + 0.614*ellRad[2]
    RAND50  <- if(accuracy) {   # take systematic location bias into account
        bias <- sqrt(sum(ctr^2)) / RAND50b
        if(bias > 2.2) {
            warning(c("\nRAND location bias estimate is ",
                      round(bias, 2), " (> 2.2),\n",
                      "more than what RAND CEP was considered for"))
        }

        ## cubic regression to take bias into account
        RAND50b * (1.0039 - 0.0528*bias + 0.4786*bias^2 - 0.0793*bias^3)
    } else {                             # ignore location bias
        RAND50b
    }                                    # if(accuracy)

    ## Grubbs-Patnaik CEP estimate based on chi^2 from Williams, 1997
    ## variance of decorrelated data (Puhek, 1992) = eigenvalues
    if(accuracy) {    # take systematic location bias into account
        GPm <- sum(ctr^2) + sum(eigVal)
        GPv <- 2*sum(eigVal^2) + 4*sum(ctr^2 * eigVal)
    } else {         # ignore location bias -> set ctr = (0,0)
        GPm <- sum(eigVal)
        GPv <- 2*sum(eigVal^2)
    }                                    # if(accuracy)
    GPn <- 2*GPm^2 / GPv                 # degrees of freedom chi-square
    Grubbs50 <- sqrt(GPv/(2*GPm) * qchisq(level, GPn))

    ## Rayleigh CEP estimate from Singh, 1992
    ## mDstCtr <- mean(getDistToCtr(xy)) # mean dist to center
    ## Rayleigh50 <- 0.9394*mDstCtr      # McMillan & McMillan
    ## unbiased sigma estimate
    sigmaHat   <- getRayParam(xy, type="Rayleigh", level=level, accuracy=accuracy)$sigma
    Rayleigh50 <- sigmaHat["sigma"] * sqrt(-2*log(1-level)) # Rayleigh quantile
    names(Rayleigh50) <- NULL

    ## CEP estimates for 90%, and 95%
    ## multiplication factors depend on aspect ratio
    ## see table 3 in McMillan & McMillan, 2008
    ## ratios <- seq(1, 4, by=0.5)       # aspect ratios considered in the table
    ## fac90  <- c(1.88, 1.93, 2.01, 2.11, 2.17, 2.28, 2.39)
    ## fac95  <- c(2.13, 2.25, 2.41, 2.44, 2.48, 2.76, 2.80)
    ## linear fit for the multiplication factors
    ## coef(lm(fac90 ~ ratios))          # R^2 = .987
    ## coef(lm(fac95 ~ ratios))          # R^2 = .939
    a <- c("90%"=1.6832142857, "95%"=1.9135714286)   # intercepts 90 and 95
    b <- c("90%"=0.1707142857, "95%"=0.2214285714)   # slopes 90 and 95

    CEP50l <- list(Rayleigh=Rayleigh50, Grubbs=Grubbs50, RAND=RAND50)
    
    ## only report the chosen estimates
    CEPl <- lapply(CEP50l[names(CEP50l) %in% type], function(x) {
                   c("50%"=x, (a + b*aspRat)*x) })
    res  <- lapply(CEPl, function(x) {
                   rbind(unit=x,
                          MOA=getMOA(x, dst=dstTarget, conversion=conversion)) })

    return(c(res, ellShape=list(c(aspectRatio=aspRat, flattening=flat)),
             ctr=list(ctr)))
}
