\name{getConfEll}
\alias{getConfEll}
\alias{getConfEll.data.frame}
\alias{getConfEll.default}
\title{Confidence ellipse}
\description{Calculates the 2D-confidence ellipse under the assumption of bivariate normality. Also includes the ellipse based on a robust estimate for the covariance matrix of the (x,y)-coordinates.}
\usage{getConfEll(xy, level = 0.5, dstTarget = 100, conversion = 'm2cm', doRob=TRUE)

\method{getConfEll}{data.frame}(xy, level = 0.5, dstTarget = 100, conversion = 'm2cm', doRob=TRUE)

\method{getConfEll}{default}(xy, level = 0.5, dstTarget = 100, conversion = 'm2cm', doRob=TRUE)}
\arguments{
  \item{xy}{(x,y)-coordinates of n points: either a numerical (n x 2)-matrix (1 row of coordinates per point), or a data frame with either the variables \code{X}, \code{Y} or \code{Point.X}, \code{Point.Y}.}
  \item{level}{a numerical value with the coverage for the confidence ellipse.}
  \item{dstTarget}{a numerical value with the distance to the target - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{doRob}{logical: should robust covariance matrix estimate be used as well?}
}
\details{No coordinate transforms are done (unlike in \code{\link[shotGroups]{groupLocation}}, \code{\link[shotGroups]{groupShape}}, \code{\link[shotGroups]{groupSpread}}), i.e., origin is not assumed to be in top-left corner, and points are not taken relative to point of aim. \cr \cr
Robust estimate for the covariance matrix of (x,y)-coordinates is from \code{\link[robustbase]{covMcd}} using the MCD algorithm. See \code{\link[shotGroups]{getCEP}} for estimates of the circular error probable.}
\value{
A list with the confidence ellipse measures.
\item{ctr}{(x,y)-coordinates group center.}
\item{ctrRob}{(x,y)-coordinates robust estimate of group center.}
\item{cov}{covariance matrix.}
\item{covRob}{robust estimate of covariance matrix.}
\item{size}{a numerical matrix with the length of the semi-major and semi-minor axis of the ellipse (in original measurement units and MOA).}
\item{sizeRob}{a numerical matrix with the length of the semi-major and semi-minor axis of the ellipse based on a robust estimate for the covariance matrix of (x,y)-coordinates (in original measurement units and MOA).}
\item{shape}{a numerical vector with the aspect ratio and the flattening of the ellipse as well as the trace and determinant of the covariance matrix.}
\item{shapeRob}{a numerical vector with the aspect ratio and the flattening of the ellipse as well as the trace and determinant based on a robust estimate for the covariance matrix of (x,y)-coordinates.}
\item{magFac}{magnification factor used to turn the error ellipse into the confidence ellipse as determined by the F(2, n-1)-distribution.}
}
\seealso{
\code{\link[shotGroups]{getCEP}},
\code{\link[robustbase]{covMcd}},
\code{\link[shotGroups]{drawEllipse}}
}
\examples{
# coordinates given by a suitable data frame
(ce <- getConfEll(DFsavage, level=0.5, dstTarget=100,
                  conversion='yd2in'))

# plot points, center, parametric confidence
# ellipse, and its axes
plot(Point.Y ~ Point.X, data=DFsavage, asp=1, pch=16)
drawEllipse(ce, axes=TRUE, fg='blue', colCtr='blue', lwd=2, pch=4, cex=2)

# coordinates given by a matrix
\dontrun{
xy <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
getConfEll(xy, level=0.5, dstTarget=25, conversion='m2cm')
}
}
