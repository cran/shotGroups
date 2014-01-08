\name{groupLocation}
\alias{groupLocation}
\alias{groupLocation.data.frame}
\alias{groupLocation.default}
\title{Accuracy: Location measures for a single group of bullet holes}
\description{Calculates location measures for a single group of bullet holes.}
\usage{groupLocation(xy, level = 0.95, plots = c('0', '1', '2'), dstTarget = 100,
              conversion = 'm2cm', target = c('BDS25m', 'DSB25m'), caliber=9)

\method{groupLocation}{data.frame}(xy, level = 0.95, plots = c('0', '1', '2'), dstTarget = 100,
              conversion = 'm2cm', target = c('BDS25m', 'DSB25m'), caliber=9)

\method{groupLocation}{default}(xy, level = 0.95, plots = c('0', '1', '2'), dstTarget = 100,
              conversion = 'm2cm', target = c('BDS25m', 'DSB25m'), caliber=9)}
\arguments{
  \item{xy}{either a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point), or a data frame with either the variables \code{X}, \code{Y} or \code{Point.X}, \code{Point.Y} as well as \code{Aim.X}, \code{Aim.Y} giving the point of aim. If missing, point of aim is assumed to be in (0,0).}
  \item{level}{a numerical value giving the level for the confidence intervals for the center (x,y)-coordinates.}
\item{plots}{a numerical value indicating whether a 2D-scatterplot should be shown. 0: no plot, 1: regular scatterplot, 2: scatterplot with added target in the background (this feature is incomplete and experimental).}
  \item{dstTarget}{a numerical value giving the distance to the target - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{target}{a character string indicating the target type to be drawn in the scatterplot. Currently, only two types of pistol targets are supported.}
  \item{caliber}{a numerical value indicating the bullet diameter in mm. Only used with \code{plots='2'}.}
}
\value{
A list with the results from the numerical analyses and statistical tests.
\item{ctr}{(x,y)-offset of group center relative to point of aim.}
\item{ctrXci}{parametric (t) and bootstrap (basic, BCa) confidence intervals for center x-coordinate (1499 replicates).}
\item{ctrYci}{parametric (t) and bootstrap (basic, BCa) confidence intervals for center y-coordinate (1499 replicates).}
\item{ctrRob}{robust estimate of group center offset relative to point of aim (MCD algorithm).}
\item{distPOA}{distance from group center to point of aim (in original measurement units and MOA).}
\item{distPOArob}{distance from robust estimate of group center to point of aim (in original measurement units and MOA).}
\item{Hotelling}{Hotelling's T^2-Test result from testing if group center equals point of aim.}
}
\seealso{
\code{\link[shotGroups]{getMOA}},
\code{\link[robustbase]{covMcd}},
\code{\link[stats]{anova.mlm}},
\code{\link[boot]{boot}},
\code{\link[boot]{boot.ci}}
}
\examples{
# coordinates given by a suitable data frame
res <- groupLocation(DFsavage, dstTarget=100, conversion='m2mm',
                     level=0.95, plots=2, target='BDS', caliber=7.62)
names(res)
res$ctr
res$distPOA
res$ctrXci
res$ctrYci

# coordinates given by a matrix
\dontrun{
# assume data from pistol shooting at 25m with 9mm ammo
# metric units
xy <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
groupLocation(xy, dstTarget=25, conversion='m2cm', plots=2,
              target='BDS25m', caliber=9)
}
}
