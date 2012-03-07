\name{groupLocation}
\alias{groupLocation}
\title{Location measures for a single group of bullet holes}
\description{Calculates location measures for a single group of bullet holes.}
\usage{groupLocation(xy, plots = c('0', '1', '2'), conversion = 'm2cm',
              dstTarget = 25, target = c('BDS', 'DSB'),
              unit=c('cm', 'in'), caliber=0.9)}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point).}
  \item{plots }{a numerical value indicating whether a 2D-scatterplot should be shown. 0: no plot, 1: regular scatterplot, 2: scatterplot with added target in the background.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{dstTarget }{a numerical value giving the distance to the target in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{target }{a character string indicating the target type to be drawn in the scatterplot. Currently, only two types of pistol targets are supported.}
  \item{unit }{a character string indicating the measurement unit for the (x,y)-coordinates. Only used with \code{plot=2}.}
  \item{caliber }{a numerical value indicating the bullet diameter in the measurement unit specified by \code{unit}. Only used with \code{plot=2}.}
}
\value{
A list with the results from the numerical analyses and statistical tests.
\item{ctr }{(x,y)-offset of group center relative to point of aim.}
\item{ctrRob }{robust estimate of group center offset relative to point of aim (MCD algorithm).}
\item{distPOA }{distance from group center to point of aim in original measurement units and MOA.}
\item{distPOArob }{distance from robust estimate of group center to point of aim in original measurement units and MOA.}
\item{Hotelling }{Hotelling's T^2-Test result from testing if group center equals point of aim.}
\item{CItX }{95\% t-confidence interval for x-coordinates.}
\item{CItY }{95\% t-confidence interval for y-coordinates.}
\item{CIbootX }{95\% bootstrap-confidence intervals for x-coordinates (percentile and BCa, 1499 replicates).}
\item{CIbootY }{95\% bootstrap-confidence intervals for y-coordinates (percentile and BCa, 1499 replicates).}
}
\seealso{
\code{\link[shotGroups]{getMOA}},
\code{\link[robustbase]{covMcd}},
\code{\link[stats]{anova.mlm}},
\code{\link[boot]{boot}},
\code{\link[boot]{boot.ci}}
}
\examples{
# assume data from pistol shooting at 25m with 9mm ammo
# metric units
xy  <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
res <- groupLocation(xy, conversion='m2cm', dstTarget=25, plots=2,
                     target='BDS', unit='cm', caliber=0.9)
names(res)
res$ctr
res$distPOA
res$CIbootX
res$CIbootY
}
