getCEP <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    UseMethod("getCEP")
}

getCEP.data.frame <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getCEP")
}

getCEP.default <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    if(!is.matrix(xy))     { stop("xy must be a matrix") }
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(!is.numeric(level)) { stop("level must be numeric") }
    if(level <= 0)         { stop("level must be > 0") }

    type <- match.arg(type, choices=c("CorrNormal", "GrubbsPearson", "GrubbsLiu",
                      "GrubbsPatnaik", "Rayleigh", "Ethridge", "RAND"), several.ok=TRUE)

    ## check if CI level is given in percent
    if(level >= 1) {
        while(level >= 1) { level <- level / 100 }
        warning(c("level must be in (0,1) and was set to ", level))
    }

    ## check if we can do robust estimation if so required
    if(nrow(xy) < 4) {
        haveRob <- FALSE
        if(doRob) {
            warning("We need >= 4 points for robust estimations")
        }
    } else {
        haveRob <- TRUE
    }                                    # if(nrow(xy) < 4)

    #####-----------------------------------------------------------------------
    ## some basic calculations used later
    if(doRob && haveRob) {        # center
        rob   <- robustbase::covMcd(xy, cor=FALSE)
        ctr   <- rob$center                       # robust estimate: group center
        sigma <- rob$cov
    } else {
        ctr   <- colMeans(xy)
        sigma <- cov(xy)
    }

    ## make sure eigenvalues >= 0 when very small
    ev     <- eigen(sigma)$values        # eigenvalues
    eigVal <- ev*sign(ev)

    ## error ellipse characteristics -> radii = sqrt of eigenvalues
    ## aspect ratio of ellipse = sqrt of kappa condition index
    aspRat <- sqrt(kappa(sigma, exact=TRUE))
    flat   <- 1 - (1/aspRat)             # flattening

    #####-----------------------------------------------------------------------
    ## CEP estimate based on correlated bivariate normal distribution -> corrNorm.R
    CorrNorm <- if(accuracy) {
        ## quantile from offset circle probability
        qmvnEll(level, mu=numeric(ncol(xy)), sigma=sigma, e=diag(ncol(xy)), x0=ctr)
    } else {
        if(ncol(xy) == 2) {              # exact Hoyt distribution
            HP <- getHoytParam(sigma)
            qHoyt(level, qpar=HP$q, omega=HP$omega)
        } else {
            qmvnEll(level, mu=numeric(ncol(xy)), sigma=sigma, e=diag(ncol(xy)),
                    x0=numeric(ncol(xy)))
        }
    }

    names(CorrNorm) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Pearson CEP estimate based on Pearson three-moment central
    ## chi^2 approximation (Grubbs, 1964, p55-56)
    ## variance of decorrelated data = eigenvalues (Puhek, 1992)
    GPP <- getGrubbsParam(sigma, ctr=ctr, accuracy=accuracy)
    GrubbsPearson <- qChisqGrubbs(level, m=GPP$m, v=GPP$v, nPrime=GPP$nPrime, type="Pearson")
    names(GrubbsPearson) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Patnaik CEP estimate based on Patnaik two-moment central
    ## chi^2 approximation (Grubbs, 1964, p54)
    GrubbsPatnaik <- qChisqGrubbs(level, m=GPP$m, v=GPP$v, n=GPP$n, type="Patnaik")
    names(GrubbsPatnaik) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Liu CEP estimate based on four-moment non-central chi^2
    ## approximation (Liu, Tang & Zhang, 2009)
    GrubbsLiu <- qChisqGrubbs(level, m=GPP$m, v=GPP$v, muX=GPP$muX, varX=GPP$varX,
                              l=GPP$l, delta=GPP$delta, type="Liu")
    names(GrubbsPatnaik) <- NULL
    
    #####-----------------------------------------------------------------------
    ## Rayleigh CEP estimate from Williams, 1997
    RayParam <- getRayParam(xy, accuracy=accuracy)
    Rayleigh <- qRayleigh(level, scale=RayParam$sigma["sigma"])
    names(Rayleigh) <- NULL
    if((aspRat > 4) && ("Rayleigh" %in% type)) {
        warning(c("Aspect ratio of error ellipse is ",
                 round(aspRat, 2) , " (> 4),\n",
                 "probably more than what Rayleigh CEP should be considered for"))
    }

    #####-----------------------------------------------------------------------
    ## Ethridge CEP estimate from Ethridge (1983) after Puhek (1992)
    lnDTC <- if(accuracy) {              # log distance to group center (radius)
        rSqSum <- sqrt(rowSums(xy^2))    # log radii to origin = point of aim
        log(rSqSum)
    } else {
        log(getDistToCtr(xy))            # log radii to group center
    }

    mLnDTC   <- mean(lnDTC)              # mean log radius
    medLnDTC <- median(lnDTC)            # median log radius
    varLnDTC <- var(lnDTC)               # variance log radius

    ## weighted mean after Hogg (1967)
    ## sample kurtosis log radius
    kLnDTC <- mean((lnDTC - mLnDTC)^4) / mean((lnDTC - mLnDTC)^2)^2
    dHogg  <- pmax(1 + (0.03 * (kLnDTC-3)^3 * (lnDTC-medLnDTC)^2 / varLnDTC), 0.01)
    wHogg  <- (1/dHogg) / sum(1/dHogg)   # weighting factors
    uHogg  <- sum(wHogg * lnDTC)         # log median radius estimate
    Ethridge50 <- exp(uHogg)
    Ethridge   <- if(level != 0.5) {
        if("Ethridge" %in% type) {
            warning("Ethridge CEP estimate is only available for level 0.5")
        }
        NA
    } else {
        Ethridge50
    }

    #####-----------------------------------------------------------------------
    ## modified RAND-234 CEP estimate for 50% from Williams, 1997
    ## using the semi-major and semi-minor axes of the error ellipse (PCA)
    RAND <- if(ncol(xy) == 2) {          # only available for 2D case
        RAND50MPI <- 0.563*sqrt(eigVal[1]) + 0.614*sqrt(eigVal[2])
        RAND50 <- if(accuracy) {        # take systematic location bias into account
            bias <- sqrt(sum(ctr^2)) / RAND50MPI
            if((bias > 2.2) && ("RAND" %in% type)) {
                warning(c("RAND location bias estimate is ",
                          round(bias, 2), " (> 2.2),\n",
                          "more than what RAND CEP should be considered for"))
            }

            ## cubic regression to take bias into account
            RAND50MPI * (1.0039 - 0.0528*bias + 0.4786*bias^2 - 0.0793*bias^3)
        } else {                         # ignore location bias
            RAND50MPI
        }                                # if(accuracy)

        if((aspRat > 4) && ("RAND" %in% type)) {
            warning(c("Aspect ratio of error ellipse is ",
                      round(aspRat, 2) , " (> 4),\n",
                      "probably more than what RAND CEP should be considered for"))
        }

        ## RAND is only available for levels 0.5, 0.9, 0.95
        if(level != 0.5) {
            if("RAND" %in% type) {
                warning("RAND CEP estimate is only available for level 0.5")
            }
            NA
        } else {
            RAND50
        }   ## else if(level %in% c(0.9, 0.95)) {
            ## RAND CEP estimates for 90% and 95% from McMillan & McMillan, 2008, table 3
            ## multiplication factors depend on aspect ratio
            ## ratios <- seq(1, 4, by=0.5)       # aspect ratios considered in the table
            ## fac90  <- c(1.88, 1.93, 2.01, 2.11, 2.17, 2.28, 2.39)
            ## fac95  <- c(2.13, 2.25, 2.41, 2.44, 2.48, 2.76, 2.80)
            ## linear fit for the multiplication factors
            ## coef(lm(fac90 ~ ratios))          # R^2 = .987
            ## coef(lm(fac95 ~ ratios))          # R^2 = .939
            ## a <- c("0.9"=1.6832142857, "0.95"=1.9135714286)   # intercepts 90 and 95
            ## b <- c("0.9"=0.1707142857, "0.95"=0.2214285714)   # slopes 90 and 95
            ## ((a + b*aspRat)*RAND50)[as.character(level)]
            ## }
    } else {
        if("RAND" %in% type) {
            warning("RAND CEP estimate is only available for 2D-data")
        }
        NA
    }

    names(RAND) <- NULL

#     #####-----------------------------------------------------------------------
#     ## Valstar CEP estimate for 50% from Williams, 1997
#     ## using the semi-major and semi-minor axes of the error ellipse (PCA)
#     Valstar <- if(ncol(xy) == 2) {       # only available for 2D case
#         ValstarMPI <- if((1/aspRat) <= 0.369) {
#             0.675*sqrt(eigVal[1]) + sqrt(eigVal[2])/(1.2*sqrt(eigVal[1]))
#         } else {
#             0.562*sqrt(eigVal[1]) + 0.615*sqrt(eigVal[2])  # almost RAND
#         }
#
#         Valstar50 <- if(accuracy) {          # take systematic location bias into account
#             Valstar <- sqrt(ValstarMPI + sum(ctr^2))
#         } else {                             # ignore location bias
#             ValstarMPI
#         }                                    # if(accuracy)
#
#         ## Valstar is only available for level 0.5
#         if(level != 0.5) {
#             if("Valstar" %in% type) {
#                 warning("Valstar CEP estimate is only available for level 0.5")
#             }
#             NA
#         } else {
#             Valstar50
#         }
#     } else {
#         if("Valstar" %in% type) {
#             warning("Valstar CEP estimate is only available for 2D-data")
#         }
#         NA
#     }
#
#     names(Valstar) <- NULL

    #####-----------------------------------------------------------------------
    ## only report the chosen estimates
    CEP <- c(CorrNormal=CorrNorm, GrubbsPearson=GrubbsPearson,
             GrubbsPatnaik=GrubbsPatnaik, GrubbsLiu=GrubbsLiu, Rayleigh=Rayleigh,
             Ethridge=Ethridge, RAND=RAND)[type]

    CEPmat <- sapply(CEP, makeMOA, dst=dstTarget, conversion=conversion)

    return(c(CEP=list(CEPmat), ellShape=list(c(aspectRatio=aspRat, flattening=flat)),
             ctr=list(ctr)))
}
