## return target definition including ring radii
getTarget <-
function(x, unit="cm", dstTarget=100, conversion="m2cm") {
    UseMethod("getTarget")
}

getTarget.character <-
function(x, unit="cm", dstTarget=100, conversion="m2cm") {
    if(!(x %in% names(targets))) {
        stop("Target unknown, see help(targets) for a list")
    }

    x <- targets[[x]]
    NextMethod("getTarget")
}

getTarget.default <-
function(x, unit="cm", dstTarget=100, conversion="m2cm") {
    unit <- match.arg(unit, choices=c("cm", "mm", "m", "in", "ft", "yd", "MOA", "SMOA", "milrad"))

    ## add ring radii
    x$ringR <- with(x, c(ringD10i/2, ringD10/2,  # ring 10 with inner sub-division
                         (ringD10/2) + seq(ringW, ringW*(nRings-1), by=ringW)))

    ## infer distance unit from conversion
    unitDst <- getUnits(conversion, first=TRUE)  # unit for distance

    ## convert ring width and radii
    x$unitConv <- unit                           # add unit converted to
    if(unit %in% c("MOA", "SMOA", "milrad"))  {  # angular size
        ringConv    <- with(x, paste0(unitDst, "2", unitTarget))
        x$ringD10u  <- with(x, getMOA(ringD10,  conversion=ringConv, dst=dstTarget, type=unit))
        x$ringD10iu <- with(x, getMOA(ringD10i, conversion=ringConv, dst=dstTarget, type=unit))
        x$ringRu    <- with(x, getMOA(ringR,    conversion=ringConv, dst=dstTarget, type=unit))
        x$ringWu    <- with(x, getMOA(ringW,    conversion=ringConv, dst=dstTarget, type=unit))

        if(!is.null(x$extra)) {
            x$extraU <- with(x, getMOA(extra, conversion=ringConv, dst=dstTarget, type=unit))
        }
    } else {                                     # absolute size
        ringConv    <- with(x, getConvFac(paste0(unitTarget, "2", unit)))
        x$ringD10u  <- with(x, ringConv*ringD10)
        x$ringD10iu <- with(x, ringConv*ringD10i)
        x$ringRu    <- with(x, ringConv*ringR)
        x$ringWu    <- with(x, ringConv*ringW)

        if(!is.null(x$extra)) {
            x$extraU <- with(x, ringConv*extra)
        }
    }

    x                                    # return selected target
}

## draw a target
drawTarget <-
function(x, unit="cm", dstTarget=100, conversion="m2cm", add=FALSE, cex=par("cex")) {
    unit <- match.arg(unit, choices=c("cm", "mm", "m", "in", "ft", "yd", "MOA", "SMOA", "milrad"))

    ## get chosen target including measures converted to unit
    target <- getTarget(x, unit=unit, dstTarget=dstTarget, conversion=conversion)

    ## open plot or add to existing plot?
    if(!add) {
        xLims <- with(target, range(c(-ringRu, ringRu)))
        yLims <- with(target, range(c(-ringRu, ringRu)))
        plot(0, 0, type="n", xlim=xLims, ylim=yLims, xlab=NA, ylab=NA, asp=1)
    }                                    # if(add)

    ## do we need a special drawing function?
    if(!is.null(target$draw)) {
        fun <- eval(target$draw)
        fun(target, cex=cex)
    } else {
        drawTarget_default(target, cex=cex)
    }

    return(invisible(target))
}

## draw a target
drawTarget_default <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x, symbols(rep(0, length(ringRu)), rep(0, length(ringRu)), add=TRUE,
                    circles=rev(ringRu), bg=rev(cols), fg=rev(colsTxt), inches=FALSE))
     # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers except for bullseye (ring number maxVal)
    ## bullseye has inner sub-division -> start numbers on ring 3
    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal-1, length.out=nRings-1)) # left side of center
    rings2 <- c(rings1, rev(rings1))                         # both sides
    pos1   <- with(x, ringRu[3:length(ringRu)] - (ringWu/2)) # right side of center
    pos2   <- c(-rev(pos1), pos1)                            # both sides
    cols1  <- with(x, colsTxt[3:length(ringRu)])             # right side of center
    cols2  <- c(rev(cols1), cols1)                           # both sides

    text(pos2, 0, cex=cex, label=rings2, col=cols2)
    text(0, pos2, cex=cex, label=rings2, col=cols2)
}

#####-----------------------------------------------------------------------
## ISSF special targets
drawTarget_ISSF25mRF <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x, symbols(rep(0, length(ringRu)), rep(0, length(ringRu)), add=TRUE,
                    circles=rev(ringRu), bg=rev(cols), fg=rev(colsTxt), inches=FALSE))
    with(x, rect(-ringRu[length(ringRu)],  -extraU/2, -ringRu[4]+(0.15*ringWu), extraU/2, col="white", border=NA))
    with(x, rect( ringRu[4]-(0.15*ringWu), -extraU/2,  ringRu[length(ringRu)],  extraU/2, col="white", border=NA))
    # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers except for bullseye (ring number maxVal)
    ## bullseye has inner sub-division -> start numbers on ring 3
    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal-1, length.out=nRings-1)) # left side of center
    rings2 <- c(rings1, rev(rings1))                         # both sides
    pos1   <- with(x, ringRu[3:length(ringRu)] - (ringWu/2)) # right side of center
    pos2   <- c(-rev(pos1), pos1)                            # both sides
    cols1  <- with(x, colsTxt[3:length(ringRu)])             # right side of center
    cols2  <- c(rev(cols1), cols1)                           # both sides

    text(0, pos2, cex=cex, label=rings2, col=cols2)
}

drawTarget_ISSF50m <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x, symbols(rep(0, length(ringRu)), rep(0, length(ringRu)), add=TRUE,
                    circles=rev(ringRu), bg=rev(cols), fg=rev(colsTxt), inches=FALSE))
    with(x, symbols(0, 0, add=TRUE, circles=extraU, bg=cols[1], fg=NA, inches=FALSE))
    with(x, symbols(rep(0, length(ringRu)-3), rep(0, length(ringRu)-3), add=TRUE,
                    circles=rev(ringRu)[-(1:3)], bg=rev(cols)[-(1:3)], fg=rev(colsTxt)[-(1:3)], inches=FALSE))

    # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers 1-8
    ## bullseye has inner sub-division -> start numbers on ring 4
    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal-2, length.out=nRings-2)) # left side of center
    rings2 <- c(rings1, rev(rings1))                         # both sides
    pos1   <- with(x, ringRu[4:length(ringRu)] - (ringWu/2)) # right side of center
    pos2   <- c(-rev(pos1), pos1)                            # both sides
    cols1  <- with(x, colsTxt[4:length(ringRu)])             # right side of center
    cols2  <- c(rev(cols1), cols1)                           # both sides

    text(pos2, 0, cex=cex, label=rings2, col=cols2)
    text(0, pos2, cex=cex, label=rings2, col=cols2)
}

## draw ISSF 300m target
drawTarget_ISSF300m <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x, symbols(rep(0, length(ringRu)), rep(0, length(ringRu)), add=TRUE,
                    circles=rev(ringRu), bg=rev(cols), fg=rev(colsTxt), inches=FALSE))
    # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers except for bullseye (ring number maxVal)
    ## bullseye has inner sub-division -> start numbers on ring 3
    rings <- with(x, seq(from=maxVal-nRings+1, to=maxVal-1, length.out=nRings-1)) # left side of center
    cols  <- with(x, colsTxt[length(ringRu):3])              # right side of center
    pos   <- with(x, -ringRu[length(ringRu):3] + (ringWu/2)) # right side of center
    ang   <- x$txtRot
    angX  <- with(x, txtRot*pi/180)
    posLU <- cbind(pos, 0) %*% with(x, cbind(c(cos(  angX), sin(  angX)), c(-sin(  angX), cos(  angX))))
    posRU <- cbind(pos, 0) %*% with(x, cbind(c(cos(3*angX), sin(3*angX)), c(-sin(3*angX), cos(3*angX))))
    posRL <- cbind(pos, 0) %*% with(x, cbind(c(cos(5*angX), sin(5*angX)), c(-sin(5*angX), cos(5*angX))))
    posLL <- cbind(pos, 0) %*% with(x, cbind(c(cos(7*angX), sin(7*angX)), c(-sin(7*angX), cos(7*angX))))

    text(posLU[ , 1], posLU[ , 2], cex=cex, label=rings, col=cols, srt=  ang)
    text(posRU[ , 1], posRU[ , 2], cex=cex, label=rings, col=cols, srt=7*ang)
    text(posRL[ , 1], posRL[ , 2], cex=cex, label=rings, col=cols, srt=5*ang)
    text(posLL[ , 1], posLL[ , 2], cex=cex, label=rings, col=cols, srt=3*ang)
}

#####-----------------------------------------------------------------------
## NRA HPR targets
drawTarget_NRA_HPR <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x, symbols(rep(0, length(ringRu)), rep(0, length(ringRu)), add=TRUE,
         circles=rev(ringRu), bg=rev(cols), fg=rev(colsTxt), inches=FALSE))
    # abline(v=0, h=0, col="lightgray")    # add point of aim

    ## add ring numbers
    ## bullseye has inner sub-division -> start numbers on ring 2
    rings1 <- with(x, seq(from=maxVal-nRings+1, to=maxVal, length.out=nRings)) # left side of center
    outRng <- 1:(length(rings1)-2)
    outPos <- 3:length(rings1)
    rings2 <- c(rings1[outRng], rev(rings1[outRng]))         # both sides
    pos1   <- with(x, ringRu[2:length(ringRu)] - (ringWu/2)) # right side of center
    pos2   <- c(-rev(pos1[outPos]), pos1[outPos])            # both sides
    cols1  <- with(x, colsTxt[2:length(ringRu)])             # right side of center
    cols2  <- c(rev(cols1[outPos]), cols1[outPos])           # both sides

    text(pos2, 0,      cex=cex, label=rings2,      col=cols2)
    text(0, pos1[1:2], cex=cex, label=rev(rings1)[1:2], col=cols1[1:2])
    text(0, 0, cex=cex, label="X", col=cols1[1])
}

#####-----------------------------------------------------------------------
## BDMP special targets
drawTarget_BDMP3 <-
function(x, cex=par("cex")) {
    ## draw circles first, start from the outside
    with(x, symbols(rep(0, length(ringRu)), rep(0, length(ringRu)), add=TRUE,
         circles=rev(ringRu), bg=rev(cols), fg=rev(colsTxt), inches=FALSE))
    with(x, rect(-extraU/2, -extraU/2, extraU/2, extraU/2, col="black", border=NA))
    with(x, symbols(rep(0, length(ringRu)), rep(0, length(ringRu)), add=TRUE,
         circles=rev(ringRu), bg=NA, fg=rev(colsTxt), inches=FALSE))
    with(x, symbols(0, 0, add=TRUE, circles=ringRu[1], bg=cols[1], fg=NA, inches=FALSE))
    # abline(v=0, h=0, col="lightgray")    # add point of aim
}

drawTarget_BDMP4 <-
function(x, cex=par("cex")) {
    ## background black square
    with(x, rect(-extraU/2, -extraU/2, extraU/2, extraU/2, col="black", border=NA))
    ## draw circles first, start from the outside
    with(x, symbols(rep(0, length(ringRu)), rep(0, length(ringRu)), add=TRUE,
         circles=rev(ringRu), bg=rev(cols), fg=rev(colsTxt), inches=FALSE))
    # abline(v=0, h=0, col="lightgray")    # add point of aim
}
