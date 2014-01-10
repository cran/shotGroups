\name{groupSpread}
\alias{groupSpread}
\alias{groupSpread.data.frame}
\alias{groupSpread.default}
\title{Precision: Spread measures of a single group of bullet holes}
\description{Provides spread measures and their graphical representations for a single group of bullet holes.}
\usage{groupSpread(xy, plots = TRUE, level = 0.95,
            CEPtype = c('Rayleigh', 'Grubbs', 'RAND'),
            sigmaType=c('Rayleigh', 'Gauss'),
            dstTarget = 100, conversion = 'm2cm')

\method{groupSpread}{data.frame}(xy, plots = TRUE, level = 0.95,
            CEPtype = c('Rayleigh', 'Grubbs', 'RAND'),
            sigmaType=c('Rayleigh', 'Gauss'),
            dstTarget = 100, conversion = 'm2cm')

\method{groupSpread}{default}(xy, plots = TRUE, level = 0.95,
            CEPtype = c('Rayleigh', 'Grubbs', 'RAND'),
            sigmaType=c('Rayleigh', 'Gauss'),
            dstTarget = 100, conversion = 'm2cm')}
\arguments{
  \item{xy}{either a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point), or a data frame with either the variables \code{X}, \code{Y} or \code{Point.X}, \code{Point.Y} as well as \code{Aim.X}, \code{Aim.Y} giving the point of aim. If missing, point of aim is assumed to be in (0,0).}
  \item{plots}{a logical value indicating whether diagrams should be shown.}
  \item{level}{a numerical value giving the level for the confidence intervals (for standard deviations, sigma, RSD, MR) as well as for the coverage of the confidence ellipse.}
  \item{CEPtype}{string indicating which CEP estimate to report from \code{\link[shotGroups]{getCEP}}.}
  \item{sigmaType}{string indicating which estimate to use for precision parameter sigma used in radial standard deviation RSD and mean radius MR.}
  \item{dstTarget}{a numerical value giving the distance to the target - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
}
\details{
Explanations and formula for many reported precision measures like CEP, sigma, RSD, MR, FoM can be found in the references.\cr
Robust estimate for the covariance matrix of (x,y)-coordinates is from \code{\link[robustbase]{covMcd}} using the MCD algorithm. \cr
In addition to the numerical results listed below, this function produces the following diagrams:
\itemize{
\item a scatterplot of the (x,y)-coordinates together with group center, circle with average distance to center, and 100*\code{level}\%-confidence ellipse - the latter also based on a robust estimate for the covariance matrix
\item a scatterplot of the (x,y)-coordinates together with the bounding box, minimum bounding box, minimum enclosing circle, and maximum group spread
\item a histogram for distances to group center including a fitted Rayleigh distribution as well as a nonparametric kernel density estimate
}
}
\value{
A list with the results from the numerical analyses and statistical tests.
\item{sdXY}{standard deviations of x- and y coordinates (in original measurement units and MOA).}
\item{sdXci}{parametric (chi^2) and bootstrap (basic, BCa) confidence intervals for the standard deviation of x-coordinates (1499 replicates, in original measurement units and MOA).}
\item{sdYci}{parametric (chi^2) and bootstrap (basic, BCa) confidence intervals for the standard deviation of y-coordinates (1499 replicates, in original measurement units and MOA).}
\item{sdXYrob}{robust standard deviations of x- and y-coordinates (in original measurement units and MOA).}
\item{covXY}{covariance matrix of (x,y)-coordinates.}
\item{covXYrob}{robust estimate of covariance matrix of (x,y)-coordinates.}
\item{distToCtr}{mean and median distance from points to their center as well as estimated Rayleigh parameters sigma (precision), radial standard deviation RSD, and mean radius MR (in original measurement units and MOA).}
\item{sigmaCI}{parametric (chi^2) and bootstrap (basic, BCa) confidence intervals for sigma (1499 replicates, in original measurement units and MOA).}
\item{RSDci}{parametric (chi^2) and bootstrap (basic, BCa) confidence intervals for radial standard deviation RSD (1499 replicates, in original measurement units and MOA).}
\item{MRci}{parametric (chi^2) and bootstrap (basic, BCa) confidence intervals for mean radius MR (1499 replicates, in original measurement units and MOA).}
\item{maxPairDist}{maximum pairwise distance between points (center-to-center, a.k.a. maximum spread, in original measurement units and MOA).}
\item{groupRect}{width and height of bounding box with diagonal and figure of merit FoM (average side length, in original measurement units and MOA).}
\item{groupRectMin}{width and height of minimum-area bounding box with diagonal and figure of merit FoM (average side length, in original measurement units and MOA).}
\item{minCircleRad}{radius for the minimum enclosing circle (in original measurement units and MOA).}
\item{confEll}{length of semi-major and semi-minor axis of the confidence ellipse (in original measurement units and MOA).}
\item{confEllRob}{length of semi-major and semi-minor axis of the confidence ellipse based on a robust estimate for the covariance matrix (in original measurement units and MOA).}
\item{confEllShape}{aspect ratio and flattening of the confidence ellipse.}
\item{confEllShapeRob}{aspect ratio and flattening of the confidence ellipse based on a robust estimate for the covariance matrix.}
\item{CEP}{estimate(s) for the circular error probable (CEP, in original measurement units and MOA).}
}
\references{
\url{http://ballistipedia.com/index.php?title=Describing_Precision}

\url{http://ballistipedia.com/index.php?title=Measuring_Precision}
}
\seealso{
\code{\link[shotGroups]{getDistToCtr}},
\code{\link[shotGroups]{getMaxPairDist}},
\code{\link[shotGroups]{getBoundingBox}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getMinCircle}},
\code{\link[shotGroups]{getConfEll}},
\code{\link[shotGroups]{getCEP}},
\code{\link[shotGroups]{getRayParam}},
\code{\link[shotGroups]{getMOA}},
\code{\link[graphics]{hist}},
\code{\link[boot]{boot}},
\code{\link[boot]{boot.ci}},
\code{\link[stats]{kernel}},
\code{\link[robustbase]{covMcd}},
}
\examples{
# coordinates given by a suitable data frame
res <- groupSpread(DFtalon, CEPtype=c("Rayleigh", "Grubbs"), level=0.95,
                   sigmaType='Rayleigh', dstTarget=10, conversion='m2mm')
names(res)
res$sdXYrob
res$distToCtr
res$maxPairDist
res$CEP

# coordinates given by a matrix
\dontrun{
xy <- matrix(round(rnorm(200, 0, 5), 2), ncol=2)
groupSpread(xy, level=0.5, dstTarget=25, conversion='m2cm')
}
}
