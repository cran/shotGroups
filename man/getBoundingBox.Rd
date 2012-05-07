\name{getBoundingBox}
\alias{getBoundingBox}
\title{Bounding box for a set of 2D-points}
\description{Calculates the vertices of the (axis-parallel) bounding box given a set of 2D-coordinates.}
\usage{getBoundingBox(xy)}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point).}
}
\value{
A list containing the rectangle's vertices, width, and height.
  \item{pts }{a numerical 4-vector giving the coordinates \code{xleft}, \code{ybottom}, \code{xright}, \code{ytop} that can be passed to \code{\link[graphics]{rect}}.}
  \item{width }{width of the box.}
  \item{height }{height of the box.}
}
\seealso{
\code{\link[shotGroups]{drawBox}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getMinCircle}}
}
\examples{
xy <- matrix(round(rnorm(20, 100, 15), 1), ncol=2)
bb <- getBoundingBox(xy)
plot(xy, pch=16)
drawBox(bb$pts[1], bb$pts[2], bb$pts[3], bb$pts[4],
        fg='blue', colCtr='blue', pch=4, cex=2)
}
