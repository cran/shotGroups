#####---------------------------------------------------------------------------
## multivariate normal offset ellipse/circle probabilities
#####---------------------------------------------------------------------------

## e characterizes the integration ellipsoid: (x-x0)' e (x-x0) < r^2 -> e = S^{-1}
## cdf
pmvnEll <-
function(r=1, sigma=diag(2), mu, e, x0, lower.tail=TRUE) {
    if(missing(mu)) { mu <- numeric(ncol(sigma)) }
    if(missing(x0)) { x0 <- numeric(ncol(sigma)) }
    if(missing(e))  { e  <- diag(ncol(sigma)) }

    ## check e, sigma positive definite
    eEV <- eigen(e)$values
    sEV <- eigen(sigma)$values
    if(!all(eEV >= -sqrt(.Machine$double.eps) * abs(eEV[1]))) {
        stop("e is numerically not positive definite")
    }

    if(!all(sEV >= -sqrt(.Machine$double.eps) * abs(sEV[1]))) {
        stop("sigma is numerically not positive definite")
    }

    ## check dimensions match
    stopifnot(length(x0) == length(mu),
              length(x0) == ncol(e),
              length(x0) == ncol(sigma))

    pp   <- numeric(length(r))               # initialize probabilities to 0
    keep <- which((r >= 0) & is.finite(r))   # keep non-negative x
    if(length(keep) < 1) {                   # nothing to do
        pp[!is.finite(r)] <- NA
        return(pp)
    }

    ## 1: Mahalanobis transform wrt e -> integrate over unit disc
    L      <- chol(e, pivot=TRUE)
    Esqrt  <- L[ , order(attr(L, "pivot"))]
    xmu1   <- Esqrt %*% (x0-mu)
    sigma1 <- Esqrt %*% sigma %*% t(Esqrt)

    ## 2: rotate with eigenvectors of sigma1 -> decorrelate normal
    S1eig <- eigen(sigma1)
    xmu2  <- t(S1eig$vectors) %*% xmu1

    ## non-centrality parameters
    ncp <- xmu2^2 / S1eig$values
    cqf <- sapply(r[keep], function(x) {
        CompQuadForm::farebrother(x^2, lambda=S1eig$values, delta=ncp)$res })

    ## CompQuadForm returns survival probabilities (1-F)
    pp[keep] <- if(lower.tail) { 1-cqf } else { cqf }

    ## special cases not caught so far
    pp[!is.finite(r)]    <- NA
    pp[which(r == -Inf)] <- 0
    pp[which(r ==  Inf)] <- 1

    return(pp)
}

## quantile function
qmvnEll <-
function(p, sigma=diag(2), mu, e, x0, lower.tail=TRUE, loUp=NULL) {
    if(missing(mu)) { mu <- numeric(ncol(sigma)) }
    if(missing(x0)) { x0 <- numeric(ncol(sigma)) }
    if(missing(e))  { e  <- diag(ncol(sigma)) }

    ## check e, sigma positive definite
    eEV <- eigen(e)$values
    sEV <- eigen(sigma)$values
    if(!all(eEV >= -sqrt(.Machine$double.eps) * abs(eEV[1]))) {
        stop("e is numerically not positive definite")
    }

    if(!all(sEV >= -sqrt(.Machine$double.eps) * abs(sEV[1]))) {
        stop("sigma is numerically not positive definite")
    }

    ## check dimensions match
    stopifnot(length(x0) == length(mu),
              length(x0) == ncol(e),
              length(x0) == ncol(sigma))

    qq   <- as.numeric(rep(NA, length(p)))
    keep <- which((p >= 0) & (p < 1))
    if(length(keep) < 1) { return(qq) }

    ## determine search interval(s) for uniroot()
    if(is.null(loUp)) {                  # no search interval given
        ## use Grubbs-Liu chi^2 quantile +- 300% for root finding
        GP <- getGrubbsParam(sigma=sigma, ctr=(x0-mu), accuracy=TRUE)
        qGrubbs <- qChisqGrubbs(p[keep], m=GP$m, v=GP$v, n=GP$n, muX=GP$muX,
                                varX=GP$varX, l=GP$l, delta=GP$delta,
                                lower.tail=lower.tail, type="Liu")
        qLo  <- max(0, qGrubbs - 3*qGrubbs)
        qUp  <- qGrubbs + 3*qGrubbs
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

    cdf <- function(r, p, x0, e, mu, sigma, lower.tail) {
        pmvnEll(r=r, sigma=sigma, mu=mu, e=e, x0=x0, lower.tail=lower.tail) - p
    }

    getQ <- function(p, x0, e, mu, sigma, loUp, lower.tail) {
        tryCatch(uniroot(cdf, interval=loUp, p=p, x0=x0, e=e, mu=mu,
                         sigma=sigma, lower.tail=lower.tail)$root,
                 error=function(e) return(NA))
    }

    qq[keep] <- unlist(Map(getQ, p=p[keep], x0=list(x0), e=list(e),
                           mu=list(mu), sigma=list(sigma), loUp=loUp[keep],
						   lower.tail=lower.tail[1]))

    return(qq)
}
