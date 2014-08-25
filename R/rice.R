#####---------------------------------------------------------------------------
## Rice distribution
## http://reference.wolfram.com/language/ref/RiceDistribution.html
## uncorrelated bivariate normal distribution with equal variances
## rewritten in polar coordinates
## pdf, cdf, and inverse cdf of the distribution of the radius around an
## offset center: nu is the offset (distance to origin), sigma the scale
#####---------------------------------------------------------------------------

## estimate Rice parameters nu, sigma, MR, RSD, from set of 2D coordinates
getRiceParam <-
function(xy, level=0.95, doRob=FALSE, type=c("LiZhangDai", "MOM")) {
    UseMethod("getRiceParam")
}

getRiceParam.data.frame <-
function(xy, level=0.95, doRob=FALSE, type=c("LiZhangDai", "MOM")) {
    xy <- getXYmat(xy)
    NextMethod("getRiceParam")
}

getRiceParam.default <-
function(xy, level=0.95, doRob=FALSE, type=c("LiZhangDai", "MOM")) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("x must be (n x 2)-matrix") }

    type <- match.arg(type)

    ## check if we can do robust estimation if so required
    if(nrow(xy) < 4) {
        haveRob <- FALSE
        if(doRob) {
            warning("We need >= 4 points for robust estimations")
        }
    } else {
        haveRob <- TRUE
        rob     <- robustbase::covMcd(xy, cor=FALSE)
    }                                    # if(nrow(xy) < 4)

    ctr <- if(doRob && haveRob) {        # estimated center
        rob$center                       # robust estimate
    } else {
        colMeans(xy)                     # coord-wise mean
    }                                    # if(doRob && haveRob)

    ## get sigma estimate from Rayleigh distribution
    sigmaHat <- getRayParam(xy=xy, level=level, doRob=doRob)$sigma

    ## estimate nu^2 -> but E(xBar'xBar) = mu'mu + (2/N)*sigma^2
    N    <- nrow(xy)
    bias <- (2/N)*sigmaHat["sigma"]^2

    ## conventional estimator: xBar'xBar - (2/N)*sigmaHat^2
    ce <- sum(ctr^2) - bias

    nuSqHat <- if(type == "MOM") {
        ## set xBar'xBar - (2/N)*sigmaHat^2 to 0 when (2/N)*sigmaHat^2 > xBar'xBar
        ifelse(ce >= 0, ce, 0)
    } else if(type == "LiZhangDai") {
        ## ML -> bad for low SNR -> instead: Li, Zhang & Dai, 2009
        max(ce, (1/(bias+1)) * sum(ctr^2))
    }

    ## c4 correction for negative bias due to taking square root (concave)
    nuHat <- (1/c4(2*N+1))*sqrt(nuSqHat)

    ## radial mean and sd
    MSD <- getMSDfromRice(nu=nuHat, sigma=sigmaHat["sigma"])

    return(list(nu=setNames(nuHat, NULL),
             sigma=sigmaHat,
                MR=setNames(MSD$mean, NULL),
               RSD=setNames(MSD$sd,   NULL)))
}

## Laguerre half polynomial
LaguerreHalf <-
function(x) {
    a   <- -x/2
    bI0 <- exp(log(besselI(a, nu=0, expon.scaled=TRUE)) + a)
    bI1 <- exp(log(besselI(a, nu=1, expon.scaled=TRUE)) + a)
    exp(x/2) * ((1-x)*bI0 - x*bI1)
}

doubleFactorial <-
function(x, log=FALSE) {
    y   <- (x + 1)/2
    lDF <- lgamma(2*y) - (lgamma(y) + (y-1) * log(2))
    if(log) {
        lDF
    } else {
        exp(lDF)
    }
}

## mean radial error and sd of radial error from Rice parameters nu, sigma
getMSDfromRice <-
function(nu, sigma) {
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    argL  <- recycle(nu, sigma)
    nu    <- argL[[1]]
    sigma <- argL[[2]]

    L05   <- LaguerreHalf(-0.5 * nu^2 / sigma^2)
    rMean <- sigma * sqrt(pi/2) * L05
    rVar  <- 2*sigma^2 + nu^2 - (pi * sigma^2 / 2) * L05^2

    ## for large signal-to-noise ratios, use approximation (Foi, 2011)
    s2nr <- nu/sigma                           # signal to noise ratio
    rMean[which(s2nr > 52)] <- nu + sigma^2/(2*nu)
    rVar[ which(s2nr > 52)] <- sigma^2 - sigma^4/(2*nu^2)

    return(list(mean=rMean, sd=sqrt(rVar)))
}

#####---------------------------------------------------------------------------
## pdf of Rice distribution
dRice <-
function(x, nu, sigma) {
    is.na(x)     <- is.nan(x)                  # replace NaN with NA
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    argL  <- recycle(x, nu, sigma)
    x     <- argL[[1]]
    nu    <- argL[[2]]
    sigma <- argL[[3]]

    dens <- numeric(length(x))                 # initialize density to 0
    keep <- which((x >= 0) | !is.finite(x))    # keep non-negative x, NA, -Inf, Inf
    if(length(keep) < 1) { return(dens) }      # nothing to do

    lfac1 <- log(x[keep]) - 2*log(sigma[keep])
    lfac2 <-   -(x[keep]^2 + nu[keep]^2) / (2*sigma[keep]^2)
    bArg  <- abs(x[keep]   * nu[keep]   /     sigma[keep]^2)
    lfac3 <- log(besselI(bArg, nu=0, expon.scaled=TRUE)) + bArg
    res   <- exp(lfac1+lfac2+lfac3)            # this may be NaN
    dens[keep] <- ifelse(is.nan(res), 0, res)  # if so, set to 0

    ## for large signal-to-noise ratios, use normal approximation
    s2nr <- nu/sigma                           # signal to noise ratio 
    rMSD <- getMSDfromRice(nu, sigma)          # M, SD of radial error

    keepS2NR24 <- keep[keep %in% which(s2nr > 24)]
    keepS2NR52 <- keep[keep %in% which(s2nr > 52)]
    dens[keepS2NR24] <- dnorm(x[keep], mean=rMSD$mean[keep], sd=rMSD$sd[keep])
    dens[keepS2NR52] <- dnorm(x[keep], mean=nu[keep],        sd=sigma[keep])

    return(dens)
}

#####---------------------------------------------------------------------------
## cdf Rice distribution
pRice <-
function(q, nu, sigma, lower.tail=TRUE) {
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    argL  <- recycle(q, nu, sigma)
    q     <- argL[[1]]
    nu    <- argL[[2]]
    sigma <- argL[[3]]

    pp   <- numeric(length(q))               # initialize probabilities to 0
    keep <- which((q >= 0) | !is.finite(q))  # keep non-negative q, NA, NaN, -Inf, Inf

    aQ  <- nu/sigma
    b0Q <- 0
    b1Q <- q/sigma

    if(lower.tail) {
        pp[keep] <-      marcumQ(aQ[keep], b0Q[keep], nu=1) - marcumQ(aQ[keep], b1Q[keep], nu=1)
        ## special cases not caught so far
        pp[which(q == -Inf)] <- 0
        pp[which(q ==  Inf)] <- 1
    } else {
        pp[keep] <- 1 - (marcumQ(aQ[keep], b0Q[keep], nu=1) - marcumQ(aQ[keep], b1Q[keep], nu=1))
        ## special cases not caught so far
        pp[which(q < 0)]    <- 1
        pp[which(q == Inf)] <- 0
    }

    ## for large signal-to-noise ratios, use normal approximation
    s2nr <- nu/sigma                         # signal to noise ratio 
    rMSD <- getMSDfromRice(nu, sigma)        # M, SD of radial error

    keepS2NR24 <- keep[keep %in% which(s2nr > 24)]
    keepS2NR52 <- keep[keep %in% which(s2nr > 52)]
    pp[keepS2NR24] <- pnorm(q[keep], mean=rMSD$mean[keep], sd=rMSD$sd[keep], lower.tail=lower.tail)
    pp[keepS2NR52] <- pnorm(q[keep], mean=nu[keep],        sd=sigma[keep],   lower.tail=lower.tail)

    return(pp)
}

#####---------------------------------------------------------------------------
## Rice quantile function through root finding of cdf
qRice <-
function(p, nu, sigma, lower.tail=TRUE, loUp=NULL) {
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    argL  <- recycle(p, nu, sigma)
    p     <- argL[[1]]
    nu    <- argL[[2]]
    sigma <- argL[[3]]

    qq   <- as.numeric(rep(NA, length(p)))
    keep <- which((p >= 0) & (p < 1))
    if(length(keep) < 1) { return(qq) }

    if(is.null(loUp)) {                  # no search interval given
        ## use Grubbs chi^2 quantile for root finding
        ## Grubbs-Liu chi^2 and Rice can diverge
        GP <- getGPfromRP(nu, sigma)     # Grubbs parameters
        qGrubbs   <- qChisqGrubbs(p[keep], m=GP$m, v=GP$v, muX=GP$muX,
                                  varX=GP$varX, l=GP$l, delta=GP$delta,
                                  lower.tail=lower.tail, type="Liu")
        qGrubbs.6 <- qChisqGrubbs(0.6, m=GP$m, v=GP$v, muX=GP$muX,
                                  varX=GP$varX, l=GP$l, delta=GP$delta,
                                  lower.tail=lower.tail, type="Liu")
        qLo  <- ifelse(p[keep] <= 0.5, 0,         0.25*qGrubbs)
        qUp  <- ifelse(p[keep] <= 0.5, qGrubbs.6, 1.75*qGrubbs)
        loUp <- split(cbind(qLo, qUp), seq_along(p))
    } else {
        if(is.matrix(loUp)) {
            loUp <- split(loUp, 1:nrow(loUp))
        } else if(is.vector(loUp)) {
            loUp <- list(loUp)
        } else if(!is.list(loUp)) {
            stop("loUp must be a list, a matrix, a vector, or missing entirely")
        }
    }

    cdf <- function(x, p, nu, sigma, lower.tail) {
        pRice(x, nu=nu, sigma=sigma, lower.tail=lower.tail) - p
    }

    getQ <- function(p, nu, sigma, loUp, lower.tail) {
        tryCatch(uniroot(cdf, interval=loUp, p=p, nu=nu, sigma=sigma,
                         lower.tail=lower.tail)$root,
                 error=function(e) return(NA))
    }

    qq[keep] <- unlist(Map(getQ, p=p[keep], nu=nu[keep], sigma=sigma[keep],
                       loUp=loUp[keep], lower.tail=lower.tail[1]))

    ## for large signal-to-noise ratios, use normal approximation
    s2nr <- nu/sigma                         # signal to noise ratio 
    rMSD <- getMSDfromRice(nu, sigma)        # M, SD of radial error

    keepS2NR24 <- keep[keep %in% which(s2nr > 24)]
    keepS2NR52 <- keep[keep %in% which(s2nr > 52)]
    qq[keepS2NR24] <- qnorm(p[keep], mean=rMSD$mean[keep], sd=rMSD$sd[keep], lower.tail=lower.tail)
    qq[keepS2NR52] <- qnorm(p[keep], mean=nu[keep],        sd=sigma[keep],   lower.tail=lower.tail)

    return(qq)
}

#####---------------------------------------------------------------------------
## random numbers from Rice distribution
rRice <-
function(n, nu, sigma, method=c("eigen", "chol", "cdf"), loUp=NULL) {
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)
    method <- match.arg(method)

    ## if n is a vector, its length determines number of random variates
    n     <- if(length(n) > 1) { length(n) } else { n }
    nu    <- nu[1]                      # only first shape parameter is used
    sigma <- sigma[1]                   # only first scale parameter is used

    rn <- if(method == "eigen") {
        ## simulated 2D normal vectors with mean 0
        X  <- matrix(rnorm(n*2), nrow=n)    # with identity cov-mat
        xy <- X %*% diag(rep(sigma, 2))     # with cov-mat according to sigma

        ## move to mean nu
        xyMove <- sweep(xy, 2, nu, FUN="+")
        sqrt(rowSums(xyMove^2))             # distances to center
    } else if(method == "chol") {
        covMat <- diag(rep(sigma^2, 2))
        CF     <- chol(covMat, pivot=TRUE)  # Cholesky-factor
        idx    <- order(attr(CF, "pivot"))
        CFord  <- CF[, idx]

        ## simulated 2D normal vectors with mean 0 and cov-mat according to sigma
        xy <- matrix(rnorm(n*2), nrow=n) %*% CFord

        ## move to mean nu
        xyMove <- sweep(xy, 2, nu, FUN="+")
        sqrt(rowSums(xyMove^2))          # distances to center
    } else if(method == "cdf") {
        ## root finding of pRice() given uniform random probabilities:
        ## find x such that F(x) - U = 0
        cdf <- function(x, u, nu, sigma) {
            pRice(x, nu=nu, sigma=sigma) - u
        }

        ## find quantile via uniroot() with error handling
        getQ <- function(u, nu, sigma, loUp) {
            tryCatch(uniroot(cdf, interval=loUp, u=u, nu=nu, sigma=sigma)$root,
                     error=function(e) return(NA))
        }

        u <- runif(n)                        # uniform random numbers

        ## determine search interval(s) for uniroot()
        if(is.null(loUp)) {                  # no search interval given
            ## use Grubbs chi^2 quantile for root finding
            ## Grubbs-Liu chi^2 and Rice can diverge
            GP <- getGPfromRP(nu, sigma)   # Grubbs parameters and quantiles
            qGrubbs   <- qChisqGrubbs(u, m=GP$m, v=GP$v, muX=GP$muX,
                                      varX=GP$varX, l=GP$l, delta=GP$delta, type="Liu")
            qGrubbs.6 <- qChisqGrubbs(0.6, m=GP$m, v=GP$v, muX=GP$muX,
                                      varX=GP$varX, l=GP$l, delta=GP$delta, type="Liu")
            qLo  <- ifelse(u <= 0.5, 0,         0.25*qGrubbs)
            qUp  <- ifelse(u <= 0.5, qGrubbs.6, 1.75*qGrubbs)
            loUp <- split(cbind(qLo, qUp), seq_along(u))
        } else {
            if(is.matrix(loUp)) {
                loUp <- split(loUp, 1:nrow(loUp))
            } else if(is.vector(loUp)) {
                loUp <- list(loUp)
            } else if(!is.list(loUp)) {
                stop("loUp must be a list, a matrix, a vector, or missing entirely")
            }
        }

        unlist(Map(getQ, u=u, nu=nu, sigma=sigma, loUp=loUp))
    }

    return(rn)
}
