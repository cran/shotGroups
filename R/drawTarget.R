drawTarget <-
function(target=c("BDS", "DSB"), unit=c("cm", "in"),
         lty=par("lty"), lwd=par("lwd"), cex=par("cex")) {
    target <- match.arg(target)
    unit   <- match.arg(unit)

    ## conversion factor to inch
    convFac <- ifelse(unit == "in", 1/2.54, 1)

    ## ring radii in cm
    rWidth  <- convFac * 2.5
    nRings  <- 10
    ringRad <- switch(target,
        BDS = c(convFac*1.25, seq(rWidth, rWidth*nRings, by=rWidth)),
        DSB = c(convFac*1.25, seq(rWidth, rWidth*nRings, by=rWidth)))

    ## set up color pattern
    colsBDS <- rep(c(396, 24, 396), c(3, 2, 6))
    colsDSB <- rep(c(24, 396),      c(5,    6))
    cols    <- switch(target,
        BDS = colors()[colsBDS],
        DSB = colors()[colsDSB])

    ## inverted colors for rings and text
    colsBDStxt <- rep(c(24, 396, 24), c(3, 2, 6))
    colsDSBtxt <- rep(c(396, 24),     c(5,    6))
    colsTxt    <- switch(target,
        BDS = colors()[colsBDStxt],
        DSB = colors()[colsDSBtxt])

    ## draw circles first, start from the outside
    n <- length(ringRad)
    symbols(rep(0, n), rep(0, n), add=TRUE, circles=rev(ringRad),
            bg=rev(cols), fg=rev(colsTxt), inches=FALSE)
    abline(v=0, h=0, col="lightgray")      # add point of aim

    ## add ring numbers except for bullseye
    nr     <- c(1:(nRings-1), (nRings-1):1)
    pos    <- c(-(rev( ringRad[3:n]) - rWidth/2), ringRad[3:n] - rWidth/2)
    colsNr <- c(  rev(colsTxt[3:n]),              colsTxt[3:n])
    text(pos, 0, cex=cex, label=nr, col=colsNr)
    text(0, pos, cex=cex, label=nr, col=colsNr)
}

# x <- -50:50
# y <- -50:50
# dev.new()
# plot(x, y, type="n", asp=1)
# drawTarget("BDS")
#
# dev.new()
# plot(x, y, type="n", asp=1)
# drawTarget("DSB")
