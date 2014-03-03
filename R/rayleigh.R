## bias correction factor for estimate of Rayleigh sigma parameter
## http://ballistipedia.com/index.php?title=Measuring_Precision#Estimating_.CF.83
getRayParam <-
function(xy, level=0.95, accuracy=FALSE) {
    UseMethod("getRayParam")
}

getRayParam.data.frame <-
function(xy, level=0.95, accuracy=FALSE) {
    xy <- getXYmat(xy)
    NextMethod("getRayParam")
}

getRayParam.default <-
function(xy, level=0.95, accuracy=FALSE) {
    if(!is.matrix(xy))     { stop("xy must be a matrix") }
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(!is.numeric(level)) { stop("level must be numeric") }
    if(level <= 0)         { stop("level must be > 0") }

    ## check if CI level is given in percent
    if(level >= 1) {
        while(level >= 1) { level <- level / 100 }
        warning(c("level must be in (0,1) and was set to ", level))
    }

    N     <- nrow(xy)
    alpha <- 1-level

    if(accuracy) {                       # Singh C2
        rSqSum <- sum(rowSums(xy^2))     # sum squared radii

        ## bias correction factor
        ## use exp(lgamma()) because gamma() will be infinite for large N
        corrFac <- sqrt(N) * exp(lgamma(N-1) - lgamma((2*N - 1)/2))
    } else {                             # Singh C1
        xyCtr   <- scale(xy, scale=FALSE, center=TRUE)  # centered data
        rSqSum  <- sum(rowSums(xyCtr^2)) # sum squared radii
        corrFac <- sqrt(N) * exp(lgamma(N)   - lgamma((2*N + 1)/2))
        ## same: http://ballistipedia.com/index.php?title=Measuring_Precision
        ## exp(log(sqrt(N/pi)) + N*log(4) + lgamma(N+1) + lgamma(N) - lgamma(2*N+1))
    }                                    # if(accuracy)

    sigBias <- sqrt((1/(2*N))*rSqSum)    # biased sigma estimate
    corrR   <- ifelse(is.finite(corrFac), corrFac, 1)
    sigHat  <- corrR * sigBias
    sigCIlo <- corrR * sqrt(rSqSum / qchisq(1-(alpha/2), 2*(N-1)))
    sigCIup <- corrR * sqrt(rSqSum / qchisq(   alpha/2,  2*(N-1)))

    RSD   <- sigHat * sqrt((4-pi)/2)     # radial standard deviation
    MR    <- sigHat * sqrt(pi/2)         # means
    RSDci <- c(sigCIlo, sigCIup) * sqrt((4-pi)/2)
    MRci  <- c(sigCIlo, sigCIup) * sqrt(pi/2)

    return(list(sigma=c(sigma=sigHat, sigCIlo=sigCIlo, sigCIup=sigCIup),
                  RSD=c(RSD=RSD, RSDciLo=RSDci[1], RSDciUp=RSDci[2]),
                   MR=c(MR=MR, MRciLo=MRci[1], MRciUp=MRci[2])))
}

#####---------------------------------------------------------------------------
## Rayleigh distribution
## pdf based on sigma
dRayleigh <-
function(x, scale=1) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)

    args  <- recycle(x, scale)
    x     <- args[[1]]
    scale <- args[[2]]

    dens <- numeric(length(x))
    keep <- which((x >= 0) | !is.finite(x))
    if(length(keep) < 1) { return(dens) }

    dens[keep] <- exp(log(x[keep]) - 0.5*(x[keep]/scale[keep])^2 - 2*log(scale[keep]))

    ## special case not caught so far
    dens[is.infinite(x)] <- 0
    return(dens)
}

## Rayleigh pdf based on sigma^2
dRayleighSq <-
function(x, scaleSq=1) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)

    args    <- recycle(x, scaleSq)
    x       <- args[[1]]
    scaleSq <- args[[2]]

    dens <- numeric(length(x))
    keep <- which((x >= 0) | !is.finite(x))
    if(length(keep) < 1) { return(dens) }

    res <- exp(log(x[keep]) - 0.5*(x[keep]^2 / scaleSq[keep]) - log(scaleSq[keep]))
    dens[keep] <- ifelse(is.nan(res), 0, res)  # if NaN, set to 0

    ## special case not caught so far
    dens[is.infinite(x)] <- 0
    return(dens)
}

## Rayleigh cdf
pRayleigh <-
function(q, scale=1) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)

    args  <- recycle(q, scale)
    q     <- args[[1]]
    scale <- args[[2]]

    pp   <- numeric(length(q))
    keep <- which((q >= 0) | !is.finite(q))
    if(length(keep) < 1) { return(pp) }

    pp[keep] <- -expm1(-0.5 * (q[keep]/scale[keep])^2)

    ## some special values not caught before
    pp[which(q == -Inf)] <- 0
    pp[which(q ==  Inf)] <- 1
    return(pp)
}

## Rayleigh quantile function
qRayleigh <-
function(p, scale=1) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)

    args  <- recycle(p, scale)
    p     <- args[[1]]
    scale <- args[[2]]

    keep <- which((p >= 0) & (p < 1))
    qq   <- as.numeric(rep(NA, length(p)))
    if(length(keep) < 1) { return(qq) }

    qq[keep] <- scale[keep] * sqrt(-2 * log1p(-p[keep]))
    return(qq)
}

## random deviates
rRayleigh <-
function (n, scale=1) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)
    rn <- scale * sqrt(-2*log(runif(n)))
    return(rn)
}
