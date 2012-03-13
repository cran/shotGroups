\name{drawEllipse}
\alias{drawEllipse}
\title{Draw a single ellipse}
\description{Adds a single ellipse to an already opened plot.}
\usage{drawEllipse(ctr, shape, radius, nv = 100, fg = par('fg'),
            bg = NA, colCtr = NA, lty = par('lty'),
            lwd = par('lwd'), pch = par('pch'), cex = par('cex'))}
\arguments{
  \item{ctr }{a numerical 2-vector giving the (x,y)-coordinates of the center.}
  \item{shape }{a numerical symmetric (2 x 2)-matrix whose eigen-structure determines the ellipse's shape.}
  \item{radius }{a numerical value giving the ellipse's magnification factor.}
  \item{nv }{number of vertices in the approximating polygon.}
  \item{fg }{color of the ellipse's rim.}
  \item{bg }{color of the ellipse's fill color. Set to \code{NA} for a fully transparent ellipse.}
  \item{colCtr }{color of the center point. Set to \code{NA} to omit.}
  \item{lty }{line type of the ellipse.}
  \item{lwd }{line width of the ellipse.}
  \item{pch }{symbol used for the center of the ellipse.}
  \item{cex }{magnification factor for symbol used for the center of the ellipse.}
}
\details{This function is mainly a wrapper for \code{\link[graphics]{polygon}}.}
\seealso{
\code{\link[graphics]{polygon}}
}
\examples{
## error ellipse for a set of points
xy     <- matrix(round(rnorm(100, 0, 8), 2), ncol=2)
ctr    <- colMeans(xy)
covMat <- cov(xy)
plot(xy, pch=16, asp=1)
drawEllipse(ctr, covMat, radius=1, fg='blue', colCtr='blue',
            lwd=2, pch=4, cex=2)
}
