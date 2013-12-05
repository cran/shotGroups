\name{getConfEll}
\alias{getConfEll}
\title{Confidence ellipse}
\description{Gives the 2D-confidence ellipse under the assumption of bivariate normality. Also includes the ellipse based on a robust estimate for the covariance matrix of the (x,y)-coordinates.}
\usage{getConfEll(xy, level = 0.5, dstTarget = 25, conversion = 'm2cm', doRob=TRUE)}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point).}
  \item{level }{a numerical value with the coverage for the confidence ellipse.}
  \item{dstTarget }{a numerical value with the distance to the target in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{doRob }{logical: should robust covariance matrix estimate also be used?}
}
\details{Robust estimate for the covariance matrix of (x,y)-coordinates is from \code{\link[robustbase]{covMcd}} using the MCD algorithm. See \code{\link[shotGroups]{getCEP}} for estimates of the circular error probable.}
\value{
A list with the confidence ellipse measures.
\item{size }{a numerical matrix with the length of the semi-major and semi-minor axis of the ellipse (in original measurement units and MOA).}
\item{sizeRob }{a numerical matrix with the length of the semi-major and semi-minor axis of the ellipse based on a robust estimate for the covariance matrix of (x,y)-coordinates (in original measurement units and MOA).}
\item{shape }{a numerical vector with the aspect ratio and the flattening of the ellipse.}
\item{shapeRob }{a numerical vector with the aspect ratio and the flattening of the ellipse based on a robust estimate for the covariance matrix of (x,y)-coordinates.}
\item{magFac }{magnification factor used to turn the error ellipse into the confidence ellipse as determined by the F(2, n-1)-distribution.}
}
\seealso{
\code{\link[shotGroups]{getCEP}},
\code{\link[robustbase]{covMcd}}
}
\examples{
xy  <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
(ce <- getConfEll(xy, level=0.5, dstTarget=25,
                  conversion='m2cm'))

# plot points, center, parametric confidence
# ellipse, and its axes
ctr    <- colMeans(xy)            # center point
eigVec <- eigen(cov(xy))$vectors  # eigenvectors
eigScl <- eigVec \%*\% diag(ce$size["unit", ])

# matrix with scaled ellipse axes
xMat <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
yMat <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])

plot(xy, asp=1, pch=16)
drawEllipse(ctr, cov(xy), radius=ce$magFac, fg="blue")
matlines(xMat, yMat, lty=1, col="blue")
points(ctr[1], ctr[2], pch=4, col="red", cex=2)
}
