\name{drawBox2}
\alias{drawBox2}
\title{Draw an oriented box}
\description{Adds an oriented box to an already opened plot.}
\usage{drawBox2(pts, fg = par('fg'), bg = NA, colCtr = NA,
         lty = par('lty'), lwd = par('lwd'), pch = par('pch'),
         cex = par('cex'))}
\arguments{
  \item{pts }{a numerical (4 x 2)-matrix giving the (x,y)-coordinates of the ordered box vertices.}
  \item{fg }{color of the box' rim.}
  \item{bg }{color of the box' fill color. Set to \code{NA} for a fully transparent box.}
  \item{colCtr }{color of the center point. Set to \code{NA} to omit.}
  \item{lty }{line type of the box.}
  \item{lwd }{line width of the box.}
  \item{pch }{symbol used for the center of the box.}
  \item{cex }{magnification factor for symbol used for the center of the box.}
}
\details{This function is mainly a wrapper for \code{\link[graphics]{polygon}}.}
\seealso{
\code{\link[shotGroups]{getMinBBox}},
\code{\link[graphics]{polygon}}
}
\examples{
xy <- matrix(round(rnorm(20, 100, 15), 1), ncol=2)
bb <- getMinBBox(xy)
plot(xy, xlim=range(c(xy[ , 1], bb$pts[ , 1])),
         ylim=range(c(xy[ , 2], bb$pts[ , 2])), asp=1, pch=16)
drawBox2(bb$pts, fg='blue', colCtr='blue', pch=4, cex=2)
}
