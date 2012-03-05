\name{groupSpread}
\alias{groupSpread}
\title{Spread measures of a single group of bullet holes}
\description{Provides spread measures and their graphical representations for a single group of bullet holes.}
\usage{groupSpread(xy, plots = TRUE,
            conversion = 'm2cm', dstTarget = 25)}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point).}
  \item{plots }{a logical value indicating whether diagrams should be shown.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{dstTarget }{a numerical value giving the distance to the target in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
}
\details{
In addition to the numerical results listed below, this function produces the following diagrams:
\itemize{
\item a scatterplot of the (x,y)-coordinates together with group center, circle with average distance to center, and 95\%-confidence ellipse - the latter also based on a robust estimate for the covariance matrix
\item a scatterplot of the (x,y)-coordinates together with the minimum bounding box, minimum enclosing circle, and maximum group spread
\item a histogram for distances to group center including a nonparametric kernel density estimate
}
}
\value{
A list with the results from the numerical analyses and statistical tests.
\item{sdXY }{standard deviations of x- and y coordinates in original measurement units and MOA.}
\item{sdXYrob }{robust standard deviations of x- and y-coordinates in original measurement units and MOA (MCD algorithm).}
\item{covXY }{covariance matrix of (x,y)-coordinates.}
\item{covXYrob }{robust estimate of covariance matrix of (x,y)-coordinates (MCD algorithm).}
\item{corXY }{correlation matrix of (x,y)-coordinates.}
\item{corXYrob }{robust estimate of correlation matrix of (x,y)-coordinates (MCD algorithm).}
\item{meanDistToCtr }{average distance from points to their center in original measurement units and MOA.}
\item{maxPairDist }{maximum pairwise distance between points in original measurement units and MOA.}
\item{groupWidth }{width of minimum bounding in original measurement units and MOA.}
\item{groupHeight }{height of minimum bounding in original measurement units and MOA.}
\item{minCircleRad }{radius for the minimum enclosing circle in original measurement units and MOA.}
}
\seealso{
\code{\link[shotGroups]{getDistToCtr}},
\code{\link[shotGroups]{getMaxPairDist}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getMinCircle}},
\code{\link[shotGroups]{getMOA}},
\code{\link[graphics]{hist}},
\code{\link[stats]{kernel}},
\code{\link[robustbase]{covMcd}},
}
\examples{
xy  <- matrix(round(rnorm(200, 0, 5), 2), ncol=2)
res <- groupSpread(xy, conversion='m2cm', dstTarget=25)
names(res)
res$sdXYrob
res$corXY
res$meanDistToCtr
res$maxPairDist
}
