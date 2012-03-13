\name{groupSpread}
\alias{groupSpread}
\title{Spread measures of a single group of bullet holes}
\description{Provides spread measures and their graphical representations for a single group of bullet holes.}
\usage{groupSpread(xy, plots = TRUE, level = 0.5,
            dstTarget = 25, conversion = 'm2cm')}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point).}
  \item{plots }{a logical value indicating whether diagrams should be shown.}
  \item{level }{a numerical value giving the coverage for the confidence ellipse.}
  \item{dstTarget }{a numerical value giving the distance to the target in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
}
\details{
Robust estimate for the covariance matrix of (x,y)-coordinates is from \code{\link[robustbase]{covMcd}} using the MCD algorithm. \cr
In addition to the numerical results listed below, this function produces the following diagrams:
\itemize{
\item a scatterplot of the (x,y)-coordinates together with group center, circle with average distance to center, and 100*\code{level}\%-confidence ellipse - the latter also based on a robust estimate for the covariance matrix
\item a scatterplot of the (x,y)-coordinates together with the minimum bounding box, minimum enclosing circle, and maximum group spread
\item a histogram for distances to group center including a nonparametric kernel density estimate
}
}
\value{
A list with the results from the numerical analyses and statistical tests.
\item{sdXY }{standard deviations of x- and y coordinates (in original measurement units and MOA).}
\item{CIsdXY }{parametric 95\%-confidence interval for the standard deviations of x- and y coordinates (in original measurement units).}
\item{sdXYrob }{robust standard deviations of x- and y-coordinates (in original measurement units and MOA).}
\item{covXY }{covariance matrix of (x,y)-coordinates.}
\item{covXYrob }{robust estimate of covariance matrix of (x,y)-coordinates.}
\item{meanDistToCtr }{average distance from points to their center (in original measurement units and MOA).}
\item{maxPairDist }{maximum pairwise distance between points (center-to-center, a.k.a. maximum spread, in original measurement units and MOA).}
\item{groupRect }{width and height of minimum bounding box (in original measurement units and MOA).}
\item{minCircleRad }{radius for the minimum enclosing circle (in original measurement units and MOA).}
\item{confEll }{length of semi-major and semi-minor axis of the confidence ellipse (in original measurement units and MOA).}
\item{confEllRob }{length of semi-major and semi-minor axis of the confidence ellipse based on a robust estimate for the covariance matrix (in original measurement units and MOA).}
\item{confEllShape }{aspect ratio and flattening of the confidence ellipse.}
\item{confEllShapeRob }{aspect ratio and flattening of the confidence ellipse based on a robust estimate for the covariance matrix.}
\item{CEPrand }{modified RAND R-234 estimate for the circular error probable (CEP, in original measurement units and MOA).}
}
\seealso{
\code{\link[shotGroups]{getDistToCtr}},
\code{\link[shotGroups]{getMaxPairDist}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getMinCircle}},
\code{\link[shotGroups]{getConfEll}},
\code{\link[shotGroups]{getCEP}},
\code{\link[shotGroups]{getMOA}},
\code{\link[graphics]{hist}},
\code{\link[stats]{kernel}},
\code{\link[robustbase]{covMcd}},
}
\examples{
xy  <- matrix(round(rnorm(200, 0, 5), 2), ncol=2)
res <- groupSpread(xy, level=0.5, dstTarget=25, conversion='m2cm')
names(res)
res$sdXYrob
res$meanDistToCtr
res$maxPairDist
res$CEPrand
}
