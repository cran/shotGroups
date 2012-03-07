drawBox <-
function(xleft, ybottom, xright, ytop, fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
   rect(xleft, ybottom, xright, ytop, col=bg, border=fg, lty=lty, lwd=lwd)
   ctr <- c(xleft + (xright-xleft) / 2, ybottom + (ytop-ybottom) / 2)
   points(ctr[1], ctr[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawBox2 <-
function(pts, fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
   polygon(pts, col=bg, border=fg, lty=lty, lwd=lwd)
   ctr <- pts[1, ] + 0.5 * (pts[3, ] - pts[1, ])
   points(ctr[1], ctr[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawCircle <-
function(x, y, radius, nv=100, fg=par("fg"), bg=NA,
         colCtr=NA, lty=par("lty"), lwd=par("lwd"),
         pch=par("pch"), cex=par("cex")) {
    if(!is.numeric(x))      { stop("x must be numeric") }
    if(!is.numeric(y))      { stop("y must be numeric") }
    if(!is.numeric(radius)) { stop("radius must be numeric") }

    angles <- seq(0, 2*pi, length.out=nv)
    circ   <- cbind(x + radius*cos(angles), y + radius*sin(angles))

    polygon(circ[-1, ], border=NA, col=bg)
    lines(circ, col=fg, lwd=lwd, lty=lty)
    points(x, y, col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawEllipse <-
function(ctr, shape, radius=1, nv=100, fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    if(!is.numeric(ctr))   { stop("ctr must be numeric") }
    if(!is.vector(ctr))    { stop("ctr must be a vector") }
    if(length(ctr) != 2)   { stop("ctr must have length two") }
    if(!is.matrix(shape))  { stop("shape must be a matrix") }
    if(!is.numeric(shape)) { stop("shape must be numeric") }
    if((nrow(shape) != 2) | (ncol(shape) != 2)) { stop("shape must be a (2 x 2)-matrix") }
    if(!isTRUE(all.equal(max(abs(shape - t(shape))), 0))) {
        stop("shape must be symmetric")
    }

    RR     <- chol(shape, pivot=TRUE)      # Cholesky-decomposition
    angles <- seq(0, 2*pi, length.out=nv)  # angles in radians
    ell    <- radius * cbind(cos(angles), sin(angles)) %*% RR      # ellipse
    ellCtr <- sweep(ell, 2, ctr, "+")      # move ellipse to center

    ## draw center and ellipse
    points(ctr[1], ctr[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)  # center
    polygon(ellCtr, border=fg, col=bg, lwd=lwd, lty=lty)           # ellipse
}
