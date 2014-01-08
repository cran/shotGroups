getColors <-
function(n) {
    hues <- seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

drawBox <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    UseMethod("drawBox")
}

drawBox.list <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    x <- x$pts
    NextMethod("drawBox")
}

drawBox.default <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    rect(x[1], x[2], x[3], x[4], col=bg, border=fg, lty=lty, lwd=lwd)
    ctr <- c(x[1] + (x[3]-x[1]) / 2, x[2] + (x[4]-x[2]) / 2)
    points(ctr[1], ctr[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawBox2 <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    UseMethod("drawBox2")
}

drawBox2.list <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    x <- x$pts
    NextMethod("drawBox2")
}

drawBox2.default <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    polygon(x, col=bg, border=fg, lty=lty, lwd=lwd)
    ctr <- x[1, ] + 0.5 * (x[3, ] - x[1, ])
    points(ctr[1], ctr[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawCircle <-
function(x, radius, nv=100, fg=par("fg"), bg=NA,
         colCtr=NA, lty=par("lty"), lwd=par("lwd"),
         pch=par("pch"), cex=par("cex")) {
    UseMethod("drawCircle")
}

drawCircle.list <-
function(x, radius, nv=100, fg=par("fg"), bg=NA,
         colCtr=NA, lty=par("lty"), lwd=par("lwd"),
         pch=par("pch"), cex=par("cex")) {
    radius <- x$rad
    x      <- x$ctr
    NextMethod("drawCircle", radius=radius)
}

drawCircle.default <-
function(x, radius, nv=100, fg=par("fg"), bg=NA,
         colCtr=NA, lty=par("lty"), lwd=par("lwd"),
         pch=par("pch"), cex=par("cex")) {
    if(!is.numeric(x))      { stop("x must be numeric") }
    if(!is.vector(x))       { stop("x must be a vector") }
    if(length(x) != 2)      { stop("x must have length 2") }
    if(!is.numeric(radius)) { stop("radius must be numeric") }

    angles <- seq(0, 2*pi, length.out=nv)
    circ   <- cbind(x[1] + radius*cos(angles), x[2] + radius*sin(angles))

    polygon(circ[-1, ], border=NA, col=bg)
    lines(circ, col=fg, lwd=lwd, lty=lty)
    points(x[1], x[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawEllipse <-
function(x, shape, radius, nv=100, axes=FALSE,
         fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    UseMethod("drawEllipse")
}

drawEllipse.list <-
function(x, shape, radius, nv=100, axes=FALSE,
        fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    if(missing(shape))  { shape  <- x$cov }
    if(missing(radius)) { radius <- x$magFac }
    x <- x$ctr
    NextMethod("drawEllipse", shape=shape, radius=radius)
}

drawEllipse.default <-
function(x, shape, radius=1, nv=100, axes=FALSE,
         fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    if(!is.numeric(x))     { stop("x must be numeric") }
    if(!is.vector(x))      { stop("x must be a vector") }
    if(length(x) != 2)     { stop("x must have length two") }
    if(!is.matrix(shape))  { stop("shape must be a matrix") }
    if(!is.numeric(shape)) { stop("shape must be numeric") }
    if(any(dim(shape) != c(2, 2))) { stop("shape must be a (2 x 2)-matrix") }
    if(!isTRUE(all.equal(max(abs(shape - t(shape))), 0, check.attributes=FALSE))) {
        stop("shape must be symmetric")
    }

    CD     <- chol(shape, pivot=TRUE)      # Cholesky-decomposition
    CDord  <- order(attr(CD, "pivot"))
    angles <- seq(0, 2*pi, length.out=nv)  # angles in radians
    ell    <- radius * cbind(cos(angles), sin(angles)) %*% CD[ , CDord]  # ellipse
    ellCtr <- sweep(ell, 2, x, "+")        # move ellipse to center
    
    ## draw center, ellipse
    points(x[1], x[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)  # center
    polygon(ellCtr, border=fg, col=bg, lwd=lwd, lty=lty)       # ellipse

    ## draw axes
    if(axes) {
        eig    <- eigen(shape)
        eigScl <- eig$vectors %*% diag(radius * sqrt(eig$values))
        
        # matrix with scaled ellipse axes
        xMat <- rbind(x[1] + eigScl[1, ], x[1] - eigScl[1, ])
        yMat <- rbind(x[2] + eigScl[2, ], x[2] - eigScl[2, ])
        
        matlines(xMat, yMat, col=fg, lwd=lwd, lty=lty)
    }
}
