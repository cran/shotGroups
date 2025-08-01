drawGroup <-
function(xy, center=FALSE,
         xyTopLeft=TRUE, bb=FALSE, bbMin=FALSE, bbDiag=FALSE, minCirc=FALSE,
         minEll=FALSE,
         maxSpread=FALSE, meanDist=FALSE, confEll=FALSE, CEP=FALSE, ringID=FALSE,
         valueID=TRUE,
         doRob=FALSE, level=0.95, scaled=TRUE, caliber=9,
         dstTarget, conversion, unit="unit", alpha=0.5, target) {
    UseMethod("drawGroup")
}

drawGroup.data.frame <-
function(xy, center=FALSE,
         xyTopLeft=TRUE, bb=FALSE, bbMin=FALSE, bbDiag=FALSE, minCirc=FALSE,
         minEll=FALSE,
         maxSpread=FALSE, meanDist=FALSE, confEll=FALSE, CEP=FALSE, ringID=FALSE,
         valueID=TRUE,
         doRob=FALSE, level=0.95, scaled=TRUE, caliber=9,
         dstTarget, conversion, unit="unit", alpha=0.5, target) {
    ## distance to target from override or from data
    if(missing(dstTarget)) {
        dstTarget <- if(hasName(xy, "distance")) {
            xy[["distance"]]
        } else {
            NA_real_
        }
    }
    
    ## determine conversion factor from data if override is not given
    if(missing(conversion)) {
        conversion <- determineConversion(xy)
    }

    ## if target unspecified but data has unique target, use that
    if(missing(target)       &&
       hasName(xy, "target") &&
       length(unique(xy[["target"]])) == 1L) {
        target <- unique(xy[["target"]])
    }
       
    xy        <- getXYmat(xy, xyTopLeft=xyTopLeft, center=center)
    xyTopLeft <- FALSE                   # swapping Y was done in getXYmat()
    center    <- FALSE                   # centering was done in getXYmat()

    drawGroup(xy=xy, center=center, xyTopLeft=xyTopLeft, bb=bb, bbMin=bbMin,
              bbDiag=bbDiag, minCirc=minCirc, minEll=minEll,
              maxSpread=maxSpread,
              meanDist=meanDist, confEll=confEll, CEP=CEP, ringID=ringID,
              valueID=valueID,
              doRob=doRob, level=level, scaled=scaled, caliber=caliber,
              dstTarget=dstTarget, conversion=conversion, unit=unit,
              alpha=alpha, target=target)
}

drawGroup.default <-
function(xy, center=FALSE,
         xyTopLeft=TRUE, bb=FALSE, bbMin=FALSE, bbDiag=FALSE, minCirc=FALSE,
         minEll=FALSE,
         maxSpread=FALSE, meanDist=FALSE, confEll=FALSE, CEP=FALSE, ringID=FALSE,
         valueID=TRUE,
         doRob=FALSE, level=0.95, scaled=TRUE, caliber=9,
         dstTarget, conversion, unit="unit", alpha=0.5, target) {
    if(!is.matrix(xy))       { stop("xy must be a matrix") }
    if(!is.numeric(xy))      { stop("xy must be numeric") }
    if(ncol(xy) != 2L)       { stop("xy must have two columns") }
    if(!is.numeric(caliber)) { stop("caliber must be numeric") }
    if(caliber <= 0)         { stop("caliber must be > 0") }
    if(!all(is.numeric(level)))    { stop("level must be numeric") }
    if(any(level <= 0))      { stop("level must be > 0") }
    if(!is.numeric(alpha))   { stop("alpha must be numeric") }
    if((alpha < 0) || (alpha > 1)) { stop("alpha must be in [0,1]") }
    if(center) {
        warning("Centering only works for data frames, ignored here")
    }
    
    unit <- match.arg(tolower(unit),
                      choices=c("unit", "m", "cm", "mm", "yd", "ft", "in",
                                "deg", "moa", "smoa", "rad", "mrad", "mil"))
    CEP  <- as.character(CEP)
    if(CEP != "FALSE") {
        ## set CEP type to given estimator or to default
        CEPtype <- if(CEP != "TRUE") { CEP } else { "CorrNormal" }
    }

    ## check if CI level is given in percent
    levelTo01 <- function(level) {
        if(level >= 1) {
            while(level >= 1) { level <- level / 100 }
            warning(c("level must be in (0,1) and was set to ", level))
        }
        level
    }
    level <- vapply(level, levelTo01, numeric(1))

    if(missing(dstTarget)) {
        dstTarget <- NA_real_
    }
    
    if(missing(conversion)) {
        conversion <- NA_character_
    }
    
    haveTarget <- if(missing(target) || is.null(target) ||
                     is.na(target)   || (tolower(target) == "none")) {
        FALSE
    } else {
        if(length(unique(target)) > 1L) {
            warning("will use only 1st target")
        }

        target <- target[1]
        TRUE
    }

    ## check if we can do robust estimation if so required
    Npts <- nrow(xy)
    haveRobustbase <- requireNamespace("robustbase", quietly=TRUE)
    haveRob <- if(haveRobustbase && (Npts >= 4L)) {
        TRUE
    } else {
        if(doRob && (Npts < 4L)) {
            warning("We need >= 4 points for robust estimations")
        }
        
        if(doRob && !haveRobustbase) {
            warning("Please install package 'robustbase' for robust estimations")
        }
        
        FALSE
    }
    
    res <- vector("list", 0)             # empty list to later collect the results

    #####-----------------------------------------------------------------------
    ## prepare data
    #####-----------------------------------------------------------------------
    ## infer (x,y)-coord units from conversion
    unitDst <- getUnits(conversion, first=TRUE)  # unit for distance
    unitXY  <- getUnits(conversion, first=FALSE) # unit for xy-coords

    ## if origin is in top-left corner -> swap Y
    ## assumes that point of aim is in (0,0)
    if(xyTopLeft) { xy[ , 2] <- -xy[ , 2] }

    ## convert coords to different unit if requested
    xyNew <- if(unit == "unit") {        # keep unit
        unitXYnew <- unitXY              # new unit = old unit
        ## convert caliber given in mm to new unit of xy-coords
        convFac <- getConvFac(paste0("mm2", unitXYnew))
        calSize <- convFac * caliber/2
        xy                               # new xy-coords = old xy-coords
    } else if(tolower(unit) %in% c("deg", "moa", "smoa", "rad", "mrad", "mil")) {
        unitXYnew <- unit
        ## convert caliber given in mm to angular size
        calSize <- getMOA(caliber/2, dst=dstTarget, conversion=paste0(unitDst, "2mm"), type=unit)

        ## new xy-coords - make positive, get angular size and add sign back
        sign(xy) * getMOA(abs(xy), dst=dstTarget, conversion=conversion, type=unit)
    } else {                                # absolute size unit
        unitXYnew <- unit
        ## convert caliber given in mm to new unit of xy-coords
        convFac <- getConvFac(paste0("mm2", unitXYnew))
        calSize <- convFac * caliber/2

        ## new xy-coords by unit conversion
        xy2xyNew <- getConvFac(paste0(unitXY, "2", unitXYnew))
        xy2xyNew * xy
    }

    calSize <- unique(calSize)
    res$xy  <- xyNew                      # save converted xy-coords

    ## extract x, y coords
    X <- xyNew[ , 1]                     # x-coords
    Y <- xyNew[ , 2]                     # y-coords

    if(all(is.na(X)) || all(is.na(Y))) {
        if(nrow(na.omit(xy)) > 0L) {
            stop("No non-missing coordinates after unit conversion\n Please supply units via 'conversion'")    
        } else {
            stop("No non-missing coordinates")
        }
        
    }

    ## to determine axis limits later, collect all results in a vector
    axisLimsX <- numeric(0)
    axisLimsY <- numeric(0)

    ## determine regular vs. robust center and confidence ellipse
    if(haveRob && doRob) {
        rob <- robustbase::covMcd(xyNew, cor=FALSE)
        ctr <- rob$center                # robust estimate: group center
    } else {
        ctr <- colMeans(xyNew)           # group center
    }

    res$ctr <- ctr                       # save group center

    if(bb) {                             # bounding box
        bBox <- getBoundingBox(xyNew)
        res$bb <- bBox

        ## for axis limits
        axisLimsX <- c(axisLimsX, bBox$pts[c(1, 3)])
        axisLimsY <- c(axisLimsY, bBox$pts[c(2, 4)])
    }

    if(bbMin) {                          # minimum-area bounding box
        bBoxMin <- getMinBBox(xyNew)
        res$bbMin <- bBoxMin

        ## for axis limits
        axisLimsX <- c(axisLimsX, bBoxMin$pts[ , 1])
        axisLimsY <- c(axisLimsY, bBoxMin$pts[ , 2])
    }

    if(bbDiag) {                         # bounding box diagonal
        ## if no bb is chosen OR if regular bb is chosen
        if(bb || !any(c(bb, bbMin))) {
            bBox <- getBoundingBox(xyNew)  # regular bounding box
            res$bbDiag <- bBox$diag
        }

        if(bbMin) {
            res$bbMinDiag <- bBoxMin$diag
        }
    }

    if(minCirc) {                        # minimum enclosing circle
        mCirc <- getMinCircle(xyNew)
        res$minCirc <- mCirc

        ## for axis limits
        axisLimsX <- c(axisLimsX, mCirc$ctr[1] + mCirc$rad,
                                  mCirc$ctr[1] - mCirc$rad)
        axisLimsY <- c(axisLimsY, mCirc$ctr[2] + mCirc$rad,
                                  mCirc$ctr[2] - mCirc$rad)
    }

    if(minEll) {                         # minimum enclosing ellipse
        mEll <- getMinEllipse(xyNew)
        res$minEll <- mEll
        
        ## for axis limits
        axisLimsX <- c(axisLimsX,
                       mEll$ctr[1] + mEll$size[1],
                       mEll$ctr[1] - mEll$size[1])
        axisLimsY <- c(axisLimsY,
                       mEll$ctr[2] + mEll$size[1],
                       mEll$ctr[2] - mEll$size[1])
    }
    
    if(maxSpread) {                      # maximum group spread
        maxPD <- getMaxPairDist(xyNew)
        res$maxPairDist <- maxPD$d
    }

    if(meanDist) {                       # mean distance to group center
        meanDstCtr <- mean(getDistToCtr(xyNew))
        res$meanDist <- meanDstCtr
    }

    if(confEll) {                        # confidence ellipse
        if(doRob && haveRob) {
            cEll <- lapply(level, function(x) {
                getConfEll(xyNew, level=x, dstTarget=dstTarget, conversion=conversion, doRob=TRUE) })

            ## overwrite size and ctr with robust estimates
            cEll <- lapply(cEll, function(x) {
                x$ctr  <- x$ctrRob
                x$size <- x$sizeRob
                x })
        } else {
            cEll <- lapply(level, function(x) {
                getConfEll(xyNew, level=x, dstTarget=dstTarget, conversion=conversion, doRob=FALSE) })
        }

        cEllCopy <- lapply(cEll, function(x) {
            x$ctrRob   <- NULL
            x$covRob   <- NULL
            x$sizeRob  <- NULL
            x$shapeRob <- NULL
            x$size <- x$size["unit", ]

            ## for axis limits
            axisLimsX <- c(axisLimsX, x$ctr[1] + x$size["semi-major"],
                                      x$ctr[1] - x$size["semi-major"])
            axisLimsY <- c(axisLimsY, x$ctr[2] + x$size["semi-major"],
                                      x$ctr[2] - x$size["semi-major"])
            x })                         # drop RAD, MOA, SMOA, rad, mrad, mil
        res$confEll <- cEllCopy
    }

    if(CEP != "FALSE") {                 # circular error probable
        CEPres <- getCEP(xyNew, CEPlevel=level, dstTarget=dstTarget,
                         conversion=conversion, type=CEPtype, accuracy=FALSE)

        res$CEP <- vapply(CEPres$CEP, function(x) { x["unit", CEPtype] }, numeric(1) )

        ## for axis limits
        axisLimsX <- c(axisLimsX, vapply(CEPres$CEP, function(x) {
                c(CEPres$ctr[1] + x["unit", CEPtype], CEPres$ctr[1] - x["unit", CEPtype])
            }, numeric(2)))
        axisLimsY <- c(axisLimsY, vapply(CEPres$CEP, function(x) {
                c(CEPres$ctr[2] + x["unit", CEPtype], CEPres$ctr[2] - x["unit", CEPtype])
            }, numeric(2)))
    }

    #####-----------------------------------------------------------------------
    ## plot
    ## determine axis limits
    xLims <- range(c(X, axisLimsX))
    yLims <- range(c(Y, axisLimsY))

    ## set up colors
    cols1 <- c(ctr=rgb(255,   0, 255, maxColorValue=255),
                bb=rgb(228,  26,  28, maxColorValue=255),
             bbMin=rgb( 55, 126, 184, maxColorValue=255),
            bbDiag=rgb( 77, 175,  74, maxColorValue=255),
         bbMinDiag=rgb(152,  78, 163, maxColorValue=255),
           minCirc=rgb(255, 127,   0, maxColorValue=255),
            minEll=rgb(255, 162,  81, maxColorValue=255),
         maxSpread=rgb(255, 255,  51, maxColorValue=255),
          meanDist=rgb(166,  86,  40, maxColorValue=255),
           confEll=rgb(247, 129, 191, maxColorValue=255),
               CEP=rgb(153, 153, 153, maxColorValue=255))

    cols2 <- c(ctr=rgb(255,   0, 255, maxColorValue=255),
                bb=rgb(166, 206, 227, maxColorValue=255),
             bbMin=rgb( 31, 120, 180, maxColorValue=255),
            bbDiag=rgb(178, 223, 138, maxColorValue=255),
         bbMinDiag=rgb( 51, 160,  44, maxColorValue=255),
           minCirc=rgb(251, 154, 153, maxColorValue=255),
            minEll=rgb(249,  92,  92, maxColorValue=255),
         maxSpread=rgb(227,  26,  28, maxColorValue=255),
          meanDist=rgb(253, 191, 111, maxColorValue=255),
           confEll=rgb(255, 127,   0, maxColorValue=255),
               CEP=rgb(202, 178, 214, maxColorValue=255))

    ## color for bullet holes -> white with target, black otherwise
    pointCol <- if(haveTarget) {
        cols <- cols1
        trgt <- getTarget(target, unit="cm", dstTarget=dstTarget, conversion=conversion)
        if(hasName(trgt, "colPt")) {
            adjustcolor(trgt$colPt, alpha.f=alpha)
        } else {
            if(all(is.na(unlist(trgt$inUnit)))) {
                cols <- cols2
                rgb(0, 0, 0, alpha)
            } else {
                rgb(1, 1, 1, alpha)
            }
        }
    } else {
        cols <- cols2
        rgb(0, 0, 0, alpha)
    }

    ## start to collect legend information
    legText <- character(0)
    legCol  <- character(0)
    legLty  <- numeric(0)
    legLwd  <- numeric(0)
    legPch  <- numeric(0)

    ## distance to target may be heterogeneous
    dstTargetPlot <- paste(unique(round(na.omit(dstTarget))), collapse=", ")
    unitXYnew     <- unique(unitXYnew)
    unitXYnewPlot <- paste(na.omit(unitXYnew), collapse=", ")

    ## start with empty plot to set up region and axis labels
    plot(Y ~ X, asp=1, type="n", main="Group (x,y)-coordinates",
         xlim=xLims, ylim=yLims,
         sub=paste("distance:", dstTargetPlot, unitDst),
         xlab=paste0("X [", na.omit(unitXYnew), "]"),
         ylab=paste0("Y [", na.omit(unitXYnew), "]"))

    ## draw target background
    if(haveTarget) {
        res$target <- drawTarget(target, unit=unitXYnew, dstTarget=dstTarget,
                                 conversion=conversion, add=TRUE, cex=1.5)
    } else {
        abline(v=0, h=0, col="lightgray")  # add point of aim
    }

    ## draw bullet holes to scale
    if(scaled && !is.na(calSize)) {
        symbols(Y ~ X, asp=1, main="(x,y)-coordinates", add=TRUE,
                circles=rep(calSize, Npts), inches=FALSE,
                fg=rgb(0.3, 0.3, 0.3, alpha), bg=pointCol)
    } else {
        points(Y ~ X, pch=20, col=pointCol)
    }

    ## add group center and robust estimate for group center
    points(ctr[1], ctr[2], col=cols["ctr"], pch=4, lwd=2, cex=2)

    legText <- c(legText, "center")
    legCol  <- c(legCol, cols["ctr"])
    legLty  <- c(legLty, NA)
    legLwd  <- c(legLwd, 2)
    legPch  <- c(legPch, 4)

    if(bb) {                             # bounding box
        drawBox(bBox, fg=cols["bb"], lwd=2)
        legText <- c(legText, "bounding box")
        legCol  <- c(legCol, cols["bb"])
        legLty  <- c(legLty, 1)
        legLwd  <- c(legLwd, 2)
        legPch  <- c(legPch, NA)
        
        if(valueID) {
            txtPosW_X <- bBox$pts["xleft"] + 0.15*(bBox$pts["xright"] - bBox$pts["xleft"])
            txtPosW_Y <- bBox$pts["ytop"]  + strheight("1234")
            
            txtPosH_X <- bBox$pts["xleft"]   - strheight("1234")
            txtPosH_Y <- bBox$pts["ybottom"] + 0.85*(bBox$pts["ytop"]  - bBox$pts["ybottom"])
            
            text(x=txtPosW_X, y=txtPosW_Y, labels=round(bBox$width, 2),  col=cols["bb"], adj=c(0.5, 0.5))
            text(x=txtPosH_X, y=txtPosH_Y, labels=round(bBox$height, 2), col=cols["bb"], adj=c(0.5, 0.5), srt=90)
        }
    }

    if(bbMin) {                          # minimum bounding box
        drawBox2(bBoxMin, fg=cols["bbMin"], lwd=2)
        legText <- c(legText, "min bounding box")
        legCol  <- c(legCol, cols["bbMin"])
        legLty  <- c(legLty, 1)
        legLwd  <- c(legLwd, 2)
        legPch  <- c(legPch, NA)
        
        if(valueID) {
            ## angle of longer edge pointing up
            dPts   <- diff(bBoxMin$pts)
            idxMax <- which.max(rowSums(dPts^2))    # one of the longer edges
            idxMin <- which.min(rowSums(dPts^2))    # one of the shorter edges
            eMax   <- dPts[idxMax, ]
            eMin   <- dPts[idxMin, ]
            eMaxUp <- eMax * sign(eMax[2])          # rotate upwards 180 deg if necessary
            eMinUp <- eMin * sign(eMin[2])          # rotate upwards 180 deg if necessary
            degMax <- atan2(eMaxUp[2], eMaxUp[1])*180 / pi  # angle in degrees
            degMin <- atan2(eMinUp[2], eMinUp[1])*180 / pi  # angle in degrees
            
            ang     <- -bBoxMin$angle * pi/180
            rotMat  <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), nrow=2)
            bbMin90 <- t(rotMat %*% t(bBoxMin$pts))
            idxMax  <- which.max(rowSums(bbMin90^2))    # one of the longer edges
            idxMin  <- which.min(rowSums(bbMin90^2))    # one of the shorter edges

            txtPosW_X <- min(bbMin90[ , 1]) + 0.15*(max(bbMin90[ , 1]) - min(bbMin90[ , 1]))
            txtPosW_Y <- max(bbMin90[ , 2]) + strheight("1234")
            txtPosH_X <- min(bbMin90[ , 1]) - strheight("1234")
            txtPosH_Y <- min(bbMin90[ , 2]) + 0.85*(max(bbMin90[ , 2]) - min(bbMin90[ , 2]))
            
            # drawBox2(bbMin90)
            # text(x=txtPosW_X, y=txtPosW_Y, labels=round(bBoxMin$width, 2),
            #      col=cols["bbMin"], adj=c(0.5, 0.5))
            # text(x=txtPosH_X, y=txtPosH_Y, labels=round(bBoxMin$height, 2),
            #      col=cols["bbMin"], adj=c(0.5, 0.5), srt=90)
            txtPosWH_XY <- matrix(c(txtPosW_X, txtPosW_Y, txtPosH_X, txtPosH_Y), byrow=TRUE, nrow=2)
            ang         <- bBoxMin$angle * pi/180
            rotMat      <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), nrow=2)
            txtPos12_XY <- t(rotMat %*% t(txtPosWH_XY))

            text(x=txtPos12_XY[1, 1], y=txtPos12_XY[1, 2], labels=round(bBoxMin$width, 2),
                 col=cols["bbMin"], adj=c(0.5, 0.5), srt=degMax)
            text(x=txtPos12_XY[2, 1], y=txtPos12_XY[2, 2], labels=round(bBoxMin$height, 2),
                 col=cols["bbMin"], adj=c(0.5, 0.5), srt=degMin)
        }
    }

    if(bbDiag) {                         # bounding box diagonal
        ## if no bb is chosen or if regular bb is chosen
        if(bb || !any(c(bb, bbMin))) {
            segments(x0=bBox$pts["xleft"],  y0=bBox$pts["ybottom"],
                     x1=bBox$pts["xright"], y1=bBox$pts["ytop"],
                     col=cols["bbDiag"], lwd=2)
            legText <- c(legText, "bound box diag")
            legCol  <- c(legCol, cols["bbDiag"])
            legLty  <- c(legLty, 1)
            legLwd  <- c(legLwd, 2)
            legPch  <- c(legPch, NA)
            
            if(valueID) {
                dPts   <- c(bBox$pts["xright"]-bBox$pts["xleft"], bBox$pts["ytop"]-bBox$pts["ybottom"])
                dPtsU  <- dPts / sqrt(sum(dPts^2))
                dPtsUO <- c(dPtsU[2], -dPtsU[1])
                
                txtPos_X <- bBox$pts["xleft"]   + 0.75*(bBox$pts["xright"] - bBox$pts["xleft"])   + strheight("1234")*dPtsUO[1]
                txtPos_Y <- bBox$pts["ybottom"] + 0.75*(bBox$pts["ytop"]   - bBox$pts["ybottom"]) + strheight("1234")*dPtsUO[2]
                
                ## angle of longer edge pointing up
                eUp  <- dPts * sign(dPts[2])            # rotate upwards 180 deg if necessary
                deg  <- atan2(eUp[2], eUp[1])*180 / pi  # angle in degrees
                
                text(x=txtPos_X, y=txtPos_Y, labels=round(bBox$diag, 2),
                     col=cols["bbDiag"], srt=deg, adj=c(0,0))
            }
        }

        if(bbMin) {                      # minimum bounding box
            segments(x0=bBoxMin$pts[1, 1], y0=bBoxMin$pts[1, 2],
                     x1=bBoxMin$pts[3, 1], y1=bBoxMin$pts[3, 2],
                     col=cols["bbMinDiag"], lwd=2)
            legText <- c(legText, "min bound box diag")
            legCol  <- c(legCol, cols["bbMinDiag"])
            legLty  <- c(legLty, 1)
            legLwd  <- c(legLwd, 2)
            legPch  <- c(legPch, NA)

            if(valueID) {
                dPts   <- diff(bBoxMin$pts)
                dPtsU  <- dPts / sqrt(sum(dPts^2))
                dPtsUO <- c(dPtsU[2], -dPtsU[1])

                txtPos_X <- bBoxMin$pts[1, 1] + 0.75*(bBoxMin$pts[3, 1] - bBoxMin$pts[1, 1]) + strheight("1234")*dPtsUO[1]
                txtPos_Y <- bBoxMin$pts[1, 2] + 0.75*(bBoxMin$pts[3, 2] - bBoxMin$pts[1, 2]) + strheight("1234")*dPtsUO[2]
                
                ## angle of longer edge pointing up
                dPts <- diff(bBoxMin$pts[c(1, 3), ])
                eUp  <- dPts * sign(dPts[2])            # rotate upwards 180 deg if necessary
                deg  <- atan2(eUp[2], eUp[1])*180 / pi  # angle in degrees
                
                text(x=txtPos_X, y=txtPos_Y, labels=round(bBoxMin$diag, 2),
                     col=cols["bbDiag"], srt=deg, adj=c(0.5, 0.5))
            }
        }
    }

    if(minCirc) {                        # minimum enclosing circle
        drawCircle(mCirc, fg=cols["minCirc"], lwd=2)
        legText <- c(legText, "min circle")
        legCol  <- c(legCol, cols["minCirc"])
        legLty  <- c(legLty, 1)
        legLwd  <- c(legLwd, 2)
        legPch  <- c(legPch, NA)
        
        if(valueID) {
            txtPos_X_top <- 0
            txtPos_Y_top <- mCirc$rad + strheight("1234")
            ang    <- pi/4
            rotMat <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), nrow=2)
            txtPos_XY <- rotMat %*% matrix(c(txtPos_X_top, txtPos_Y_top), nrow=2)
            txtPos_X  <- mCirc$ctr[1] + txtPos_XY[1, 1]
            txtPos_Y  <- mCirc$ctr[2] + txtPos_XY[2, 1]
            text(x=txtPos_X, y=txtPos_Y, labels=round(mCirc$rad, 2),
                 srt=180*ang/pi, col=cols["minCirc"], adj=c(0.5, 0.5))
        }
    }
    
    if(minEll) {                         # minimum enclosing ellipse
        drawEllipse(mEll, fg=cols["minEll"], lwd=2)
        legText <- c(legText, "min ellipse")
        legCol  <- c(legCol, cols["minEll"])
        legLty  <- c(legLty, 1)
        legLwd  <- c(legLwd, 2)
        legPch  <- c(legPch, NA)
        
        if(valueID) {
            txtPos_X_top <- 0
            txtPos_Y_top <- mEll$size[1] + strheight("1234")
            ang    <- -mEll$shape["angle"] * (pi/180)
            rotMat <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), nrow=2)
            txtPos_XY <- rotMat %*% matrix(c(txtPos_X_top, txtPos_Y_top), nrow=2)
            txtPos_X  <- mEll$ctr[1] + txtPos_XY[1, 1]
            txtPos_Y  <- mEll$ctr[2] + txtPos_XY[2, 1]
            text(x=txtPos_X, y=txtPos_Y,
                 labels=round(mEll$size[1], 2),
                 srt=180*ang/pi, col=cols["minEll"], adj=c(0.5, 0.5))
            
            txtPos_X_top <- 0
            txtPos_Y_top <- mEll$size[2] + strheight("1234")
            ang    <- -(mEll$shape["angle"]+90) * (pi/180)
            rotMat <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), nrow=2)
            txtPos_XY <- rotMat %*% matrix(c(txtPos_X_top, txtPos_Y_top), nrow=2)
            txtPos_X  <- mEll$ctr[1] + txtPos_XY[1, 1]
            txtPos_Y  <- mEll$ctr[2] + txtPos_XY[2, 1]
            text(x=txtPos_X, y=txtPos_Y,
                 labels=round(mEll$size[2], 2),
                 srt=180*ang/pi, col=cols["minEll"], adj=c(0.5, 0.5))
        }
    }

    if(maxSpread) {                      # maximum group spread
        segments(x0=xyNew[maxPD$idx[1], 1], y0=xyNew[maxPD$idx[1], 2],
                 x1=xyNew[maxPD$idx[2], 1], y1=xyNew[maxPD$idx[2], 2],
                 col=cols["maxSpread"], lwd=2)
        legText <- c(legText, "max spread")
        legCol  <- c(legCol, cols["maxSpread"])
        legLty  <- c(legLty, 1)
        legLwd  <- c(legLwd, 2)
        legPch  <- c(legPch, NA)
        
        if(valueID) {
            m      <- xyNew[maxPD$idx, ]
            dPts   <- diff(m)
            dPtsU  <- dPts / sqrt(sum(dPts^2))
            dPtsUO <- c(dPtsU[2], -dPtsU[1])
            txtPos_X <- m[1, 1] + 0.75*dPts[1] + strheight("1234")*dPtsUO[1]
            txtPos_Y <- m[1, 2] + 0.75*dPts[2] + strheight("1234")*dPtsUO[2]

            ## angle of longer edge pointing up
            eUp  <- dPts * sign(dPts[2])            # rotate upwards 180 deg if necessary
            deg  <- atan2(eUp[2], eUp[1])*180 / pi  # angle in degrees
            
            text(x=txtPos_X, y=txtPos_Y, labels=round(maxPD$d, 2),
                 col=cols["maxSpread"], srt=deg, adj=c(0,0))
        }
    }

    if(meanDist) {                       # mean distance to center
        drawCircle(ctr, radius=meanDstCtr, fg=cols["meanDist"], lwd=2)
        legText <- c(legText, "mean dist to ctr")
        legCol  <- c(legCol, cols["meanDist"])
        legLty  <- c(legLty, 1)
        legLwd  <- c(legLwd, 2)
        legPch  <- c(legPch, NA)

        if(valueID) {
            ctr <- colMeans(xy)
            txtPos_X_top <- 0
            txtPos_Y_top <- meanDstCtr + strheight("1234")
            ang    <- -pi/4
            rotMat <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), nrow=2)
            txtPos_XY <- rotMat %*% matrix(c(txtPos_X_top, txtPos_Y_top), nrow=2)
            txtPos_X  <- ctr[1] + txtPos_XY[1, 1]
            txtPos_Y  <- ctr[2] + txtPos_XY[2, 1]

            text(x=txtPos_X, y=txtPos_Y, labels=round(meanDstCtr, 2),
                 srt=180*ang/pi, col=cols["meanDist"], adj=c(0.5, 0.5))
        }
    }

    if(confEll) {                        # confidence ellipse
        ellCtrL <- lapply(cEll, function(x) { drawEllipse(x, pch=4, fg=cols["confEll"], lwd=2) })
        legText <- c(legText, paste("conf ellipse", paste(level, collapse=" ")))
        legCol  <- c(legCol, cols["confEll"])
        legLty  <- c(legLty, 1)
        legLwd  <- c(legLwd, 2)
        legPch  <- c(legPch, NA)
        
        if(valueID) {
            idx      <- which.max(ellCtrL[[1]][ , 1])
            txtPos_X <- ellCtrL[[1]][idx, 1] + strheight("1234")
            txtPos_Y <- ellCtrL[[1]][idx, 2]
            label    <- paste(c(round(cEllCopy[[1]]$size["semi-major"], 2),
                                round(cEllCopy[[1]]$size["semi-minor"], 2)),
                              collapse="; ")
            
            text(x=txtPos_X, y=txtPos_Y, labels=label,
                 srt=-90, col=cols["confEll"], adj=c(0.5, 0.5))
        }
    }

    if(CEP != "FALSE") {                 # circular error probable
        lapply(CEPres$CEP, function(x) {
            drawCircle(CEPres$ctr, x["unit", CEPtype], fg=cols["CEP"], lwd=2) })
        
        legText <- c(legText, paste("CEP", CEPtype, paste(level, collapse=" ")))
        legCol  <- c(legCol, cols["CEP"])
        legLty  <- c(legLty, 1)
        legLwd  <- c(legLwd, 2)
        legPch  <- c(legPch, NA)
        
        if(valueID) {
            txtPos_X <- CEPres$ctr[1]
            txtPos_Y <- CEPres$ctr[2] + CEPres$CEP[[1]]["unit", ] + strheight("1234")
            text(x=txtPos_X, y=txtPos_Y, labels=signif(CEPres$CEP[[1]]["unit", ], 2),
                 col=cols["CEP"], adj=c(0.5, 0.5))
        }
    }

    if(ringID && haveTarget) {
        rc <- simRingCount(xy, target=target, caliber=caliber, unit=unique(unitXY))
        res$ringCount <- with(rc, c(count=count, max=max))
        with(rc, text(xyNew[,1], xyNew[,2], label=levels(rings)[rings],
                      adj=c(0.5, 0.5), col="darkgreen"))
    }

    ## add legend
    legend(x="bottomleft", legend=legText, col=legCol, lty=legLty,
           lwd=legLwd, pch=legPch, bg=rgb(1, 1, 1, 0.7))

    #####-----------------------------------------------------------------------
    ## invisibly return converted coords and spread indicators
    return(invisible(res))
}
