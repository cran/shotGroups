getCEP <-
function(xy, CEPlevel=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    UseMethod("getCEP")
}

getCEP.data.frame <-
function(xy, CEPlevel=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getCEP")
}

getCEP.default <-
function(xy, CEPlevel=0.5, dstTarget=100, conversion="m2cm", accuracy=FALSE,
         type="CorrNormal", doRob=FALSE) {
    if(!is.matrix(xy))        { stop("xy must be a matrix") }
    if(!is.numeric(xy))       { stop("xy must be numeric") }
    if(!is.numeric(CEPlevel)) { stop("CEPlevel must be numeric") }
    if(CEPlevel <= 0)         { stop("CEPlevel must be > 0") }

    type <- match.arg(type,
                      choices=c("CorrNormal", "GrubbsPearson", "GrubbsLiu",
                                "GrubbsPatnaik", "Rayleigh", "Krempasky",
                                "Ethridge", "RAND"), several.ok=TRUE)

    ## check if CEPlevel is given in percent
    if(CEPlevel >= 1) {
        while(CEPlevel >= 1) { CEPlevel <- CEPlevel / 100 }
        warning(c("CEPlevel must be in (0,1) and was set to ", CEPlevel))
    }

    ## check if we can do robust estimation if so required
    haveRob <- if(nrow(xy) < 4) {
        if(doRob) { warning("We need >= 4 points for robust estimations") }
        FALSE
    } else {
        TRUE
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
    lambda <- ev*sign(ev)

    ## error ellipse characteristics -> radii = sqrt of eigenvalues
    ## aspect ratio of ellipse = sqrt of condition index kappa
    aspRat <- sqrt(kappa(sigma, exact=TRUE))
    flat   <- 1 - (1/aspRat)             # flattening

    #####-----------------------------------------------------------------------
    ## CEP estimate based on correlated bivariate normal distribution
    CorrNorm <- if(accuracy) {
        ## quantile from offset circle probability -> mvnEll.R
        qmvnEll(CEPlevel, mu=numeric(ncol(xy)), sigma=sigma,
                e=diag(ncol(xy)), x0=ctr)
    } else {
        if(ncol(xy) == 2) {              # exact Hoyt distribution -> 1.R
            HP <- getHoytParam(sigma)
            qHoyt(CEPlevel, qpar=HP$q, omega=HP$omega)
        } else {                         # 1D/3D case -> mvnEll.R
            qmvnEll(CEPlevel, mu=numeric(ncol(xy)), sigma=sigma,
                    e=diag(ncol(xy)), x0=numeric(ncol(xy)))
        }
    }

    names(CorrNorm) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Pearson CEP estimate based on Pearson three-moment central
    ## chi^2 approximation (Grubbs, 1964, p55-56)
    GPP <- getGrubbsParam(sigma, ctr=ctr, accuracy=accuracy)
    GrubbsPearson <- qChisqGrubbs(CEPlevel, m=GPP$m, v=GPP$v, nPrime=GPP$nPrime, type="Pearson")
    names(GrubbsPearson) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Patnaik CEP estimate based on Patnaik two-moment central
    ## chi^2 approximation (Grubbs, 1964, p54)
    GrubbsPatnaik <- qChisqGrubbs(CEPlevel, m=GPP$m, v=GPP$v, n=GPP$n, type="Patnaik")
    names(GrubbsPatnaik) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Liu CEP estimate based on four-moment non-central chi^2
    ## approximation (Liu, Tang & Zhang, 2009)
    GrubbsLiu <- qChisqGrubbs(CEPlevel, m=GPP$m, v=GPP$v, muX=GPP$muX, varX=GPP$varX,
                              l=GPP$l, delta=GPP$delta, type="Liu")
    names(GrubbsLiu) <- NULL

    #####-----------------------------------------------------------------------
    ## Rayleigh CEP estimate from Singh, 1992
    Rayleigh <- if(ncol(xy) == 2) {      # 2D -> Rayleigh distribution
        if(accuracy) {                   # POA != POI -> Rice
            RiceParam <- getRiceParam(xy, doRob=doRob)
            qRice(CEPlevel, nu=RiceParam$nu, sigma=RiceParam$sigma["sigma"])
        } else {                         # POA = POI -> Rayleigh
            RayParam <- getRayParam(xy, doRob=doRob)
            qRayleigh(CEPlevel, scale=RayParam$sigma["sigma"])
        }
    } else if(ncol(xy) == 3) {           # 3D -> Maxwell-Boltzmann distribution
        MaxParam <- getMaxParam(xy, doRob=doRob)
        if(accuracy) {                   # offset circle probability
            ## circular covariance matrix with estimated M-B param sigma
            sigMat <- diag(rep(MaxParam$sigma["sigma"]^2, ncol(xy)))
            qmvnEll(CEPlevel, sigma=sigMat, mu=numeric(ncol(xy)),
                    x0=ctr, e=diag(ncol(xy)))
        } else {                         # no offset
            qMaxwell(CEPlevel, sigma=MaxParam$sigma["sigma"])
        }
    }

    names(Rayleigh) <- NULL
    if((aspRat > 4) && ("Rayleigh" %in% type)) {
        warning(c("Aspect ratio of error ellipse is ",
                 round(aspRat, 2) , " (> 4),\n",
                 "probably more than what Rayleigh CEP should be considered for"))
    }

    #####-----------------------------------------------------------------------
    ## Krempasky CEP estimate from Krempasky, 2003
    Krempasky50 <- if(ncol(xy) == 2) {   # 2D
        if(accuracy) {
            warning("Krempasky CEP estimate is only available for accuracy=FALSE")
            NA_real_
        } else {                         # POA = POI -> Rayleigh
            ## estimated correlation, covariance and standard deviations
            rho    <- cov2cor(sigma)[1, 2]
            covXY  <- sigma[1, 2]
            sigmaX <- sqrt(sigma[1, 1])
            sigmaY <- sqrt(sigma[2, 2])

            ## rotation angle gamma and rotation matrix
            gamma <- atan((-2*covXY + sqrt(4*covXY^2 + (sigmaY^2 - sigmaX^2)^2)) /
                          (sigmaY^2 - sigmaX^2))
            
            A <- cbind(c(cos(gamma), sin(gamma)), c(-sin(gamma), cos(gamma)))
            
            ## covariance matrix, correlation and standard deviations of rotated data
            sigmaRot  <- t(A) %*% sigma %*% A     # covariance matrix
            rhoDash   <- cov2cor(sigmaRot)[1, 2]  # correlation
            sigmaDash <- sqrt(sigmaRot[1, 1])
            stopifnot(isTRUE(all.equal(sigmaRot[1, 1], sigmaRot[2, 2]))) # supposed to be equal

            CEP00 <- sqrt(2*log(2))*sigmaDash
            C2    <- 0.5*(1 - log(2)/2)
            C4    <- C2*(log(2) - log(2)^2/4 - 0.5) + 3/8 - (9/16)*log(2) + (3/16)*log(2)^2 - (1/64)*log(2)^3 - C2^2*log(2)/2
            CEP00*(1 - 0.5*C2*rhoDash^2 - 0.5*(C4 + 0.25*C2^2)*rhoDash^4)
        }
    } else {
        if("Krempasky" %in% type) {
            warning("Krempasky CEP estimate is only available for 2D-data")
        }
        NA_real_
    }

    Krempasky <- if(CEPlevel != 0.5) {
        if("Krempasky" %in% type) {
            warning("Krempasky CEP estimate is only available for CEPlevel 0.5")
        }
        NA_real_
    } else {
        Krempasky50
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
    Ethridge   <- if(CEPlevel != 0.5) {
        if("Ethridge" %in% type) {
            warning("Ethridge CEP estimate is only available for CEPlevel 0.5")
        }
        NA_real_
    } else {
        Ethridge50
    }

    #####-----------------------------------------------------------------------
    ## modified RAND-234 CEP estimate for 50% from Williams, 1997
    ## using the semi-major and semi-minor axes of the error ellipse (PCA)
    RAND <- if(ncol(xy) == 2) {          # only available for 2D case
        RAND50MPI <- 0.563*sqrt(lambda[1]) + 0.614*sqrt(lambda[2])
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

        ## RAND is only available for CEPlevel 0.5
        if(CEPlevel != 0.5) {
            if("RAND" %in% type) {
                warning("RAND CEP estimate is only available for CEPlevel 0.5")
            }
            NA_real_
        } else {
            RAND50
        }   ## else if(CEPlevel %in% c(0.9, 0.95)) {
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
            ## ((a + b*aspRat)*RAND50)[as.character(CEPlevel)]
            ## }
    } else {
        if("RAND" %in% type) {
            warning("RAND CEP estimate is only available for 2D-data")
        }
        NA_real_
    }

    names(RAND) <- NULL

#     #####-----------------------------------------------------------------------
#     ## Valstar CEP estimate for 50% from Williams, 1997
#     ## using the semi-major and semi-minor axes of the error ellipse (PCA)
#     Valstar <- if(ncol(xy) == 2) {       # only available for 2D case
#         ValstarMPI <- if((1/aspRat) <= 0.369) {
#             0.675*sqrt(lambda[1]) + sqrt(lambda[2])/(1.2*sqrt(lambda[1]))
#         } else {
#             0.562*sqrt(lambda[1]) + 0.615*sqrt(lambda[2])  # almost RAND
#         }
#
#         Valstar50 <- if(accuracy) {          # take systematic location bias into account
#             Valstar <- sqrt(ValstarMPI^2 + sum(ctr^2))
#         } else {                             # ignore location bias
#             ValstarMPI
#         }                                    # if(accuracy)
#
#         ## Valstar is only available for CEPlevel 0.5
#         if(CEPlevel != 0.5) {
#             if("Valstar" %in% type) {
#                 warning("Valstar CEP estimate is only available for CEPlevel 0.5")
#             }
#             NA_real_
#         } else {
#             Valstar50
#         }
#     } else {
#         if("Valstar" %in% type) {
#             warning("Valstar CEP estimate is only available for 2D-data")
#         }
#         NA_real_
#     }
#
#     names(Valstar) <- NULL

    #####-----------------------------------------------------------------------
    ## only report the chosen estimates
    CEP <- c(CorrNormal=CorrNorm,
             GrubbsPearson=GrubbsPearson,
             GrubbsPatnaik=GrubbsPatnaik,
             GrubbsLiu=GrubbsLiu,
             Rayleigh=Rayleigh,
             Krempasky=Krempasky,
             Ethridge=Ethridge,
             RAND=RAND)[type]

    CEPmat <- sapply(CEP, makeMOA, dst=dstTarget, conversion=conversion)

    return(c(CEP=list(CEPmat), ellShape=list(c(aspectRatio=aspRat, flattening=flat)),
             ctr=list(ctr)))
}
