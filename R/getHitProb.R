getHitProb <-
function(xy, r=1, unit="unit", dstTarget=100, conversion="m2cm",
         accuracy=FALSE, type="CorrNormal", doRob=FALSE) {
    UseMethod("getHitProb")
}

getHitProb.data.frame <-
function(xy, r=1, unit="unit", dstTarget=100, conversion="m2cm",
         accuracy=FALSE, type="CorrNormal", doRob=FALSE) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getHitProb")
}

getHitProb.default <-
function(xy, r=1, unit="unit", dstTarget=100, conversion="m2cm",
         accuracy=FALSE, type="CorrNormal", doRob=FALSE) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(!is.numeric(r))  { stop("r must be numeric") }
    if(r <= 0)          { stop("r must be > 0") }

    unit <- match.arg(unit, choices=c("unit", "m", "cm", "mm", "yd", "ft", "in", "MOA", "SMOA", "mrad", "mil"))
    type <- match.arg(type, choices=c("CorrNormal", "GrubbsPearson", "GrubbsPatnaik", "GrubbsLiu", "Rayleigh"), several.ok=TRUE)

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
        ctr   <- rob$center              # robust estimate: group center
        sigma <- rob$cov                 # robust estimate: covariance matrix
    } else {
        ctr   <- colMeans(xy)
        sigma <- cov(xy)
    }

    ## error ellipse characteristics -> radii = sqrt of eigenvalues
    ## aspect ratio of ellipse = sqrt of kappa condition index
    aspRat <- sqrt(kappa(sigma, exact=TRUE))
    flat   <- 1 - (1/aspRat)             # flattening

    ## infer (x,y)-coord units from conversion
    unitXY <- getUnits(conversion, first=FALSE)

    ## convert r to unit of (x,y)-coordinates
    rNew <- if(unit == "unit") {         # keep unit
        r                                # new r = r
    } else if(unit %in% c("MOA", "SMOA", "mrad", "mil")) {
        fromMOA(r, dst=dstTarget, conversion=conversion, type=unit)
    } else {                             # absolute size unit
        r2rNew <- getConvFac(paste0(unit, "2", unitXY))
        r2rNew * r
    }

    #####-----------------------------------------------------------------------
    ## estimate based on correlated bivariate normal distribution
    CorrNorm <- if(accuracy) {
        ## offset circle probability -> mvnEll.R
        pmvnEll(r=rNew, sigma=sigma, mu=numeric(ncol(xy)), x0=ctr, e=diag(ncol(xy)))
    } else {
        if(ncol(xy) == 2) {              # exact Hoyt distribution -> hoyt.R
            HP <- getHoytParam(sigma)
            pHoyt(rNew, qpar=HP$q, omega=HP$omega)
        } else {                         # 1D/3D case -> mvnEll.R
            pmvnEll(r=rNew, sigma=sigma, mu=numeric(ncol(xy)), x0=ctr, e=diag(ncol(xy)))
        }
    }

    names(CorrNorm) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Pearson CEP estimate based on Pearson three-moment central
    ## chi^2 approximation (Grubbs, 1964, p55-56)
    GPP <- getGrubbsParam(sigma, ctr=ctr, accuracy=accuracy)
    GrubbsPearson <- pChisqGrubbs(rNew, m=GPP$m, v=GPP$v, nPrime=GPP$nPrime, type="Pearson")
    names(GrubbsPearson) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Patnaik CEP estimate based on Patnaik two-moment central
    ## chi^2 approximation (Grubbs, 1964, p54)
    GrubbsPatnaik <- pChisqGrubbs(rNew, m=GPP$m, v=GPP$v, n=GPP$n, type="Patnaik")
    names(GrubbsPatnaik) <- NULL

    #####-----------------------------------------------------------------------
    ## Grubbs-Liu CEP estimate based on four-moment non-central chi^2
    ## approximation (Liu, Tang & Zhang, 2009)
    GrubbsLiu <- pChisqGrubbs(rNew, m=GPP$m, v=GPP$v, muX=GPP$muX, varX=GPP$varX,
                              l=GPP$l, delta=GPP$delta, type="Liu")
    names(GrubbsLiu) <- NULL

    #####-----------------------------------------------------------------------
    ## Rayleigh estimate from Singh
    Rayleigh <- if(ncol(xy) == 2) {      # 2D -> Rayleigh/Rice distribution
        if(accuracy) {                   # POA != POI -> Rice
            RiceParam <- getRiceParam(xy, doRob=doRob)
            pRice(rNew, nu=RiceParam$nu, sigma=RiceParam$sigma["sigma"])
        } else {                         # POA = POI -> Rayleigh
            RayParam <- getRayParam(xy, doRob=doRob)
            pRayleigh(rNew, scale=RayParam$sigma["sigma"])
        }
    } else if(ncol(xy) == 3) {           # 3D -> Maxwell-Boltzmann distribution
        MaxParam <- getMaxParam(xy, doRob=doRob)
        if(accuracy) {                   # offset circle probability
            ## circular covariance matrix with estimated M-B param sigma
            sigMat <- diag(rep(MaxParam$sigma["sigma"]^2, ncol(xy)))
            pmvnEll(r=rNew, sigma=sigMat, mu=numeric(ncol(xy)), x0=ctr, e=diag(ncol(xy)))
        } else {                         # no offset
            pMaxwell(rNew, sigma=MaxParam$sigma["sigma"])
        }
    }

    names(Rayleigh) <- NULL
    if((aspRat > 4) && ("Rayleigh" %in% type)) {
        warning(c("Aspect ratio of error ellipse is ",
                  round(aspRat, 2) , " (> 4),\n",
                  "probably more than what Rayleigh distribution should be considered for"))
    }

    #####-----------------------------------------------------------------------
    ## only report the chosen estimates
    CEPinv <- c(CorrNormal=CorrNorm, GrubbsPearson=GrubbsPearson,
                GrubbsPatnaik=GrubbsPatnaik, GrubbsLiu=GrubbsLiu, Rayleigh=Rayleigh)[type]
    return(CEPinv)
}
