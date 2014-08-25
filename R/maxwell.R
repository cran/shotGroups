## Maxwell-Boltzmann distribution
## http://de.scribd.com/doc/23073705/Maxwell-Distribution-Rev3

#####---------------------------------------------------------------------------
## parameter estimation -> sigma, mean, sd
getMaxParam <-
function(xyz, level=0.95, mu, doRob=FALSE) {
    UseMethod("getMaxParam")
}

getMaxParam.data.frame <-
function(xyz, level=0.95, mu, doRob=FALSE) {
    xyz <- getXYmat(xyz)
    NextMethod("getMaxParam")
}

getMaxParam.default <-
function(xyz, level=0.95, mu, doRob=FALSE) {
    if(!is.matrix(xyz))  { stop("xy must be a matrix") }
    if(!is.numeric(xyz)) { stop("xy must be numeric") }
    if(ncol(xyz) != 3)   { stop("x must be (n x 3)-matrix") }

    ## check if CI level is given in percent
    if(level >= 1) {
        while(level >= 1) { level <- level / 100 }
        warning(c("level must be in (0,1) and was set to ", level))
    }

    ## check if we can do robust estimation if so required
    if(nrow(xyz) < 4) {
        haveRob <- FALSE
        if(doRob) {
            warning("We need >= 4 points for robust estimations")
        }
    } else {
        haveRob <- TRUE
        rob     <- robustbase::covMcd(xyz, cor=FALSE)
    }                                    # if(nrow(xyz) < 4)

    N     <- nrow(xyz)
    alpha <- 1-level

    if(!missing(mu)) {                   # true mean is known
        xyzCtr <- sweep(xyz, 2, mu, "-") # center with true mean
        rSqSum <- if(doRob && haveRob) { # sum squared radii centered data
            ## = N*trace of uncorrected covariance matrix -> N*trace((N-1)/N)*cov(xyz)
            (nrow(xyz)-1)*sum(diag(rob$cov))
        } else {
            sum(xyzCtr^2)                # sum squared radii
        }                                # if(doRob && haveRob)

        varHat  <- rSqSum/(3*N)          # unbiased variance estimate
        corrFac <- 1/c4(3*N + 1)         # c4 correction with n = 3*N + 1
        chisqDF <- 3*N                   # chi^2 degrees of freedom
    } else {                             # true center is estimated
        rSqSum <- if(doRob && haveRob) { # sum squared radii centered data
            ## = N*trace of uncorrected covariance matrix -> N*trace((N-1)/N)*cov(xyz)
            (nrow(xyz)-1)*sum(diag(rob$cov))
        } else {
            xyzCtr <- scale(xyz, scale=FALSE, center=TRUE)  # centered data
            sum(xyzCtr^2)                # sum squared radii
        }                                # if(doRob && haveRob)

        ## unbiased variance estimate including Bessel correction
        varHat  <- rSqSum/(3*(N-1))
        corrFac <- 1/c4(3*N - 2)         # c4 correction with n = 3*N - 2
        chisqDF <- 3*(N-1)               # chi^2 degrees of freedom
    }                                    # if(!missing(mu))

    sigHat  <- corrFac * sqrt(varHat)    # bias-corrected sigma estimate
    sigCIlo <- corrFac * sqrt(rSqSum / qchisq(1-(alpha/2), chisqDF))
    sigCIup <- corrFac * sqrt(rSqSum / qchisq(   alpha/2,  chisqDF))

    MR    <- sigHat * sqrt(8/pi)         # radial error mean
    RSD   <- sigHat * sqrt((3*pi-8)/pi)  # radial error standard deviation
    MRci  <- c(sigCIlo, sigCIup) * sqrt(8/pi)
    RSDci <- c(sigCIlo, sigCIup) * sqrt((3*pi-8)/pi)

    return(list(sigma=c(sigma=sigHat, sigCIlo=sigCIlo, sigCIup=sigCIup),
                 mean=c(MR=MR, MRciLo=MRci[1], MRciUp=MRci[2]),
                   sd=c(RSD=RSD, RSDciLo=RSDci[1], RSDciUp=RSDci[2])))
}

#####---------------------------------------------------------------------------
## Maxwell-Boltzmann distribution
dMaxwell <-
function(x, sigma=1) {
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    args  <- recycle(x, sigma)
    x     <- args[[1]]
    sigma <- args[[2]]

    dens <- numeric(length(x))
    keep <- which((x >= 0) | !is.finite(x))
    if(length(keep) < 1) { return(dens) }

    dens[keep] <- sqrt(2/pi) * (x^2/sigma[keep]^3)*exp(-x[keep]^2/(2*sigma[keep]^2))
    return(dens)
}

## Maxwell-Boltzmann cdf
pMaxwell <-
function(q, sigma=1, lower.tail=TRUE) {
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    args  <- recycle(q, sigma)
    q     <- args[[1]]
    sigma <- args[[2]]

    pp   <- numeric(length(q))
    keep <- which((q >= 0) | !is.finite(q))

    if(lower.tail) {
        pp[keep] <- 2*pnorm(q[keep]/sigma[keep]) - 1 -
            sqrt(2/pi) * (q[keep]/sigma[keep])*exp(-q[keep]^2/(2*sigma[keep]^2))

        ## some special values not caught before
        pp[which(q == -Inf)] <- 0
        pp[which(q ==  Inf)] <- 1
    } else {
        pp[keep] <- 1 - (2*pnorm(q[keep]/sigma[keep]) - 1 -
            sqrt(2/pi) * (q[keep]/sigma[keep])*exp(-q[keep]^2/(2*sigma[keep]^2)))

        ## some special values not caught before
        pp[which(q < 0)]    <- 1
        pp[which(q == Inf)] <- 0
    }
    
    return(pp)
}

## Maxwell-Boltzmann quantile function
qMaxwell <-
function(p, sigma=1, lower.tail=TRUE) {
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    args  <- recycle(p, sigma)
    p     <- args[[1]]
    sigma <- args[[2]]

    keep <- which((p >= 0) & (p < 1))
    qq   <- as.numeric(rep(NA, length(p)))
    if(length(keep) < 1) { return(qq) }

    # qmvnEll(p, sigma=diag(rep(sigma^2, 3)), mu=rep(0, 3), e=diag(3), x0=rep(0, 3))
    qq[keep] <- if(lower.tail) {
        sqrt(qgamma(p[keep],   shape=3/2, scale=2*sigma[keep]^2))
    } else {
        sqrt(qgamma(1-p[keep], shape=3/2, scale=2*sigma[keep]^2))
    }
    
    return(qq)
}

## random deviates
rMaxwell <-
function(n, sigma=1) {
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    ## simulate coords separately instead of matrix(rnorm(3*n), ncol=3)
    ## for correct cyclic replication of sigma
    xyz <- replicate(3, rnorm(n, mean=0, sd=sigma))
    sqrt(rowSums(xyz^2))
}
