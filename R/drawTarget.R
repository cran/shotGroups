drawTarget <-
function(x=c("BDS25m", "DSB25m"), unit=c("m", "cm", "mm", "in", "yd", "ft"),
         add=TRUE, cex=par("cex")) {
    x    <- match.arg(x)
    unit <- match.arg(unit)

    ## conversion factor from cm to required unit
    convFac <- getConvFac(paste0("cm2", unit))

    ## ring radii in cm
    rWidth  <- convFac * 2.5
    nRings  <- 10
    ringRad <- switch(x,
        BDS25m = c(convFac*1.25, seq(rWidth, rWidth*nRings, by=rWidth)),
        DSB25m = c(convFac*1.25, seq(rWidth, rWidth*nRings, by=rWidth)))

    ## set up color pattern
    colsBDS <- rep(c(396, 24, 396), c(3, 2, 6))
    colsDSB <- rep(c(24, 396),      c(5,    6))
    cols    <- switch(x,
        BDS25m = colors()[colsBDS],
        DSB25m = colors()[colsDSB])

    ## inverted colors for rings and text
    colsTxtBDS <- rep(c(24, 396, 24), c(3, 2, 6))
    colsTxtDSB <- rep(c(396, 24),     c(5,    6))
    colsTxt    <- switch(x,
        BDS25m = colors()[colsTxtBDS],
        DSB25m = colors()[colsTxtDSB])

    ## open plot or add to existing plot?
    if(!add) {
        xLims <- range(c(-ringRad, ringRad))
        yLims <- range(c(-ringRad, ringRad))
        plot(0, 0, type="n", xlim=xLims, ylim=yLims, xlab=NA, ylab=NA, asp=1)
    }                                    # if(add)

    ## draw circles first, start from the outside
    n <- length(ringRad)
    symbols(rep(0, n), rep(0, n), add=TRUE, circles=rev(ringRad),
            bg=rev(cols), fg=rev(colsTxt), inches=FALSE)
    abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers except for bullseye
    nr     <- c(1:(nRings-1), (nRings-1):1)
    pos    <- c(-(rev( ringRad[3:n]) - rWidth/2), ringRad[3:n] - rWidth/2)
    colsNr <- c(  rev(colsTxt[3:n]),              colsTxt[3:n])
    text(pos, 0, cex=cex, label=nr, col=colsNr)
    text(0, pos, cex=cex, label=nr, col=colsNr)
}
