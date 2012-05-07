\name{getMinCircle}
\alias{getMinCircle}
\title{Minimum enclosing circle for a set of 2D-points}
\description{Calculates center and radius of the minimum enclosing circle given a set of 2D-coordinates.}
\usage{getMinCircle(xy)}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n >= 2 points (1 row of coordinates per point).}
}
\details{Uses the Skyum algorithm based on the convex hull.}
\value{
A list containing the center and radius of the circle.
\item{ctr }{a numerical 2-vector giving the (x,y)-coordinates of the circle's center.}
\item{rad }{a numerical value giving the radius of the circle.}
}
\seealso{
\code{\link[shotGroups]{drawCircle}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getBoundingBox}}
}
\examples{
xy <- matrix(round(rnorm(20, 100, 15), 2), ncol=2)
mc <- getMinCircle(xy)

# determine axis limits so that circle will be visible
xLims <- mc$ctr[1] + c(-mc$rad, mc$rad)
yLims <- mc$ctr[2] + c(-mc$rad, mc$rad)
plot(xy, pch=16, asp=1, xlim=xLims, ylim=yLims)
drawCircle(mc$ctr[1], mc$ctr[2], mc$rad, fg='blue')
}
