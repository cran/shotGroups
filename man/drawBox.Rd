\name{drawBox}
\alias{drawBox}
\title{Draw an axis-aligned box}
\description{Adds an axis-aligned box to an already opened plot.}
\usage{drawBox(xleft, ybottom, xright, ytop, fg = par('fg'), bg = NA,
        colCtr = NA, lty = par('lty'), lwd = par('lwd'),
        pch = par('pch'), cex = par('cex'))}
\arguments{
  \item{xleft }{a vector (or scalar) of left x positions.}
  \item{ybottom }{a vector (or scalar) of bottom y positions.}
  \item{xright }{a vector (or scalar) of right x positions.}
  \item{ytop }{a vector (or scalar) of top y positions.}
  \item{fg }{color of the box' rim.}
  \item{bg }{color of the box' fill color. Set to \code{NA} for a fully transparent box.}
  \item{colCtr }{color of the center point. Set to \code{NA} to omit.}
  \item{lty }{line type of the box.}
  \item{lwd }{line width of the box.}
  \item{pch }{symbol used for the center of the box.}
  \item{cex }{magnification factor for symbol used for the center of the box.}
}
\details{This function is mainly a wrapper for \code{\link[graphics]{rect}}.}
\seealso{
\code{\link[shotGroups]{getBoundingBox}},
\code{\link[graphics]{rect}}
}
\examples{
xy <- matrix(round(rnorm(20, 100, 15), 1), ncol=2)
bb <- getBoundingBox(xy)
plot(xy, asp=1, pch=16)
drawBox(bb$pts[1], bb$pts[2], bb$pts[3], bb$pts[4],
        fg='blue', colCtr='blue', pch=4, cex=2)
}
