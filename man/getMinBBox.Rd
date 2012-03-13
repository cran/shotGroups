\name{getMinBBox}
\alias{getMinBBox}
\title{Minimum-area bounding box for a set of 2D-points}
\description{Calculates the vertices of the minimum-area, possibly oriented bounding box given a set of 2D-coordinates.}
\usage{getMinBBox(xy)}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n >= 2 points (1 row of coordinates per point).}
}
\value{
A list with the following information about the minimum-area bounding box:
\item{pts }{a (4 x 2)-matrix containing the coordinates of the (ordered) vertices.}
\item{width }{width of the box.}
\item{height }{height of the box.}
\item{angle }{orientation of the box' longer edge pointing up as returned by \code{\link{atan2}}, but in degree.}
}
\seealso{
\code{\link[shotGroups]{drawBox2}},
\code{\link[shotGroups]{getBoundingBox}},
\code{\link[shotGroups]{getMinCircle}}
}
\examples{
xy <- matrix(round(rnorm(16, 100, 15)), ncol=2)
bb <- getMinBBox(xy)                     # minimum bounding box
H  <- chull(xy)                          # convex hull

# plot original points, convex hull, and minimum bounding box
plot(xy, xlim=range(c(xy[ , 1], bb$pts[ , 1])),
         ylim=range(c(xy[ , 2], bb$pts[ , 2])), asp=1, pch=16)
polygon(xy[H, ], col=NA)                 # show convex hull
drawBox2(bb$pts, fg='blue', colCtr='blue', pch=4, cex=2)

bb$width * bb$height                     # box area
bb$angle                                 # box orientation
}
