## recycling for function arguments
recycle <- function(...) {
    dotArg <- list(...)
    maxLen <- max(sapply(dotArg, length))
    lapply(dotArg, rep, length=maxLen)
}

## bias correction factor for estimate of Rayleigh sigma parameter
## http://ballistipedia.com/index.php?title=Measuring_Precision#Estimating_.CF.83
getRayParam <-
function(xy, type=c("Rayleigh", "Gauss"), level=0.95, accuracy=FALSE) {
    UseMethod("getRayParam")
}

getRayParam.data.frame <-
function(xy, type=c("Rayleigh", "Gauss"), level=0.95, accuracy=FALSE) {
    xy <- getXYmat(xy)
    NextMethod("getRayParam")
}

getRayParam.default <-
function(xy, type=c("Rayleigh", "Gauss"), level=0.95, accuracy=FALSE) {
    if(!is.matrix(xy))     { stop("xy must be a matrix") }
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(ncol(xy) != 2)      { stop("xy must have two columns") }
    if(!is.numeric(level)) { stop("level must be numeric") }
    if(level <= 0)         { stop("level must be > 0") }
    
    type <- match.arg(type)
    
    ## check if CI level is given in percent
    if(level > 1) {
        while(level > 1) {
            level <- level / 100
        }
        warning(c("level must be in (0,1) and was set to ", level, "\n"))
    }

    N     <- nrow(xy)
    alpha <- 1-level

    if(type == "Rayleigh") {
        if(accuracy) {                       # Singh C2
            RsqSum <- sum(rowSums(xy^2))     # sum squared radii
            
            ## bias correction factor
            ## use exp(lgamma()) because gamma() will be infinite for large N
            corrFac <- sqrt(N) * exp(lgamma(N-1) - lgamma((2*N - 1)/2))
        } else {                             # Singh C1
            xyCtr   <- scale(xy, scale=FALSE, center=TRUE)  # centered data
            RsqSum  <- sum(rowSums(xyCtr^2)) # sum squared radii
            corrFac <- sqrt(N) * exp(lgamma(N)   - lgamma((2*N + 1)/2))
            ## same: http://ballistipedia.com/index.php?title=Measuring_Precision
            ## exp(log(sqrt(N/pi)) + N*log(4) + lgamma(N+1) + lgamma(N) - lgamma(2*N+1))
        }                                    # if(accuracy)
        
        sigBias <- sqrt(1/(2*N)*RsqSum)      # biased sigma estimate
        corrR   <- ifelse(is.finite(corrFac), corrFac, 1)
        sigHat  <- corrR * sigBias
        sigCIlo <- corrR * sqrt(RsqSum / qchisq(1-(alpha/2), 2*(N-1)))
        sigCIup <- corrR * sqrt(RsqSum / qchisq(   alpha/2,  2*(N-1)))
    } else if(type == "Gauss") {
        varHat  <- mean(diag(cov(xy)))       # variance estimate
        sigBias <- sqrt(varHat)              # biased sigma estimate
        
        ## bias correction factor
        ## use exp(lgamma()) because gamma() will be infinite for large N
        Ncorr   <- 2*(N-1)                   # effective N
        corrFac <- 1 / (sqrt(2/(Ncorr-1)) * exp(lgamma(Ncorr/2) - lgamma((Ncorr-1)/2)))
        
        corrG   <- ifelse(is.finite(corrFac), corrFac, 1)
        sigHat  <- corrG * sigBias
        sigCIlo <- corrG * sqrt(2*(N-1)*varHat / qchisq(1-(alpha/2), 2*(N-1)))
        sigCIup <- corrG * sqrt(2*(N-1)*varHat / qchisq(   alpha/2,  2*(N-1)))
    }                                    # if(type == "Rayleigh")

    RSD   <- sigHat * sqrt((4-pi)/2)     # radial standard deviation
    MR    <- sigHat * sqrt(pi/2)         # means
    RSDci <- c(sigCIlo, sigCIup) * sqrt((4-pi)/2)
    MRci  <- c(sigCIlo, sigCIup) * sqrt(pi/2)
    
    return(list(sigma=c(sigma=sigHat, sigCIlo=sigCIlo, sigCIup=sigCIup),
                  RSD=c(RSD=RSD, RSDciLo=RSDci[1], RSDciUp=RSDci[2]),
                   MR=c(MR=MR, MRciLo=MRci[1], MRciUp=MRci[2])))
}

## density of the Rayleigh distribution based on sigma
dRayleigh <-
function(x, scale=1) {
    ## bring x and scale to same length
    args  <- recycle(x, scale)
    x     <- args[[1]]
    scale <- args[[2]]
    
    keep <- x > 0                      # density = 0 for x <= 0
    logDens <- log(numeric(length(x))) # initialize results vector to log(0) = -Inf
    logDens[keep] <- log(x[keep]) - 0.5*(x[keep]/scale[keep])^2 - 2*log(scale[keep])
    exp(logDens)
}

## density of the Rayleigh distribution based on sigma^2
dRayleighSq <-
function(x, scaleSq=1) {
    ## bring x and scale to same length
    args    <- recycle(x, scaleSq)
    x       <- args[[1]]
    scaleSq <- args[[2]]
        
    keep <- x > 0                      # density = 0 for x <= 0
    logDens <- log(numeric(length(x))) # initialize results vector to log(0) = -Inf
    logDens[keep] <- log(x[keep]) - 0.5*(x[keep]^2 / scaleSq[keep]) - log(scaleSq[keep])
    exp(logDens)
}

qRayleigh <-
function(p, scale=1) {
    if (any(p <= 0) || any(p >= 1)) {
        stop("argument 'p' must be between 0 and 1")
    }

    res <- scale * sqrt(-2 * log(-p))
    res[scale <= 0] <- NaN
    res
}
