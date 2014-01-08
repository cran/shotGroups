\name{analyzeGroup}
\alias{analyzeGroup}
\title{Analysis for a single group of bullet holes}
\description{Performs a comprehensive numerical and graphical analysis of a single group of bullet holes.}
\usage{analyzeGroup(DF, xyTopLeft = TRUE, conversion = 'm2cm', bandW = 0.5,
             CEPtype = c('Rayleigh', 'Grubbs', 'RAND'))}
\arguments{
  \item{DF}{a data frame containing (at least) either the variables \code{Point.X}, \code{Point.Y} or \code{X}, \code{Y} defining the bullet holes. Variables \code{Distance} (distance to target), \code{Aim.X}, \code{Aim.Y} (point of aim) are useful - if they are missing, a warning is given and a default assumed.}
  \item{xyTopLeft}{a logical value indicating whether the origin of the absolute coordinate system is in the top-left corner. This is the default for data exported by OnTarget PC/TDS.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{bandW}{for argument \code{bandwith} of \code{\link[graphics]{smoothScatter}}.}
  \item{CEPtype}{string vector indicating which CEP estimate to report in \code{\link[shotGroups]{getCEP}}.}
}
\details{
Robust estimates for the group center and the covariance matrix of (x,y)-coordinates are from \code{\link[robustbase]{covMcd}} using the MCD algorithm. \cr
This function is a wrapper for \code{\link[shotGroups]{groupShape}}, \code{\link[shotGroups]{groupLocation}}, and \code{\link[shotGroups]{groupSpread}}. \cr
If the data is missing information about the point of aim, (0,0) is assumed. If distance to target is missing, 100 is assumed.\cr \cr
In addition to the numerical results listed below, this function produces the following diagrams:
\itemize{
\item a combined plot for multivariate outlier identification as produced by \code{\link[mvoutlier]{aq.plot}}
\item a scatterplot of the (x,y)-coordinates together with group center, circle with average distance to center, 50\%-confidence ellipse - the latter also based on a robust estimate for the covariance matrix
\item a scatterplot of the (x,y)-coordinates together with the minimum bounding box, minimum enclosing circle, and maximum group spread
\item a chi-square Q-Q-plot for eyeballing multivariate normality as produced by \code{\link[mvoutlier]{chisq.plot}}, including a reference line with intercept 0 and slope 1
\item a heatmap of a 2D-kernel density estimate for the (x,y)-coordinates as produced by \code{\link[graphics]{smoothScatter}} together with group center and error ellipse based on a robust estimate for the covariance matrix
\item a Q-Q-plot for x-coordinates for eyeballing normality
\item a Q-Q-plot for y-coordinates for eyeballing normality
\item a histogram for x-coordinates including a fitted normal distribution as well as a nonparametric kernel density estimate
\item a histogram for y-coordinates including a fitted normal distribution as well as a nonparametric kernel density estimate
\item a histogram for distances to group center including a fitted Rayleigh distribution as well as a nonparametric kernel density estimate
}
}
\value{
A list with the results from the numerical analyses and statistical tests.
\item{corXY}{correlation matrix of (x,y)-coordinates.}
\item{corXYrob}{robust estimate of correlation matrix of (x,y)-coordinates.}
\item{Outliers}{a vector of row indices for observations identified as outliers.}
\item{ShapiroX}{Shapiro-Wilk-Test result for normality of x-coordinates.}
\item{ShapiroY}{Shapiro-Wilk-Test result for normality of y-coordinates.}
\item{multNorm}{E-statistic-Test result for multivariate normality of (x,y)-coordinates.}

\item{sdXY}{standard deviations of x- and y coordinates (in original measurement units and MOA).}
\item{sdXci}{parametric and bootstrap confidence intervals for the standard deviation of x-coordinates (1499 replicates, in original measurement units and MOA).}
\item{sdYci}{parametric and bootstrap confidence intervals for the standard deviation of y-coordinates (1499 replicates, in original measurement units and MOA).}
\item{sdXYrob}{robust standard deviations of x- and y-coordinates (in original measurement units and MOA).}
\item{covXY}{covariance matrix of (x,y)-coordinates.}
\item{covXYrob}{robust estimate of covariance matrix of (x,y)-coordinates.}
\item{distToCtr}{mean and median distance from points to their center as well as estimated Rayleigh parameters sigma (precision), radial standard deviation RSD, and mean radius MR (in original measurement units and MOA).}
\item{sigmaCI}{parametric and bootstrap confidence intervals for sigma (1499 replicates, in original measurement units and MOA).}
\item{RSDci}{parametric and bootstrap confidence intervals for radial standard deviation RSD (1499 replicates, in original measurement units and MOA).}
\item{MRci}{parametric and bootstrap confidence intervals for mean radius MR (1499 replicates, in original measurement units and MOA).}
\item{maxPairDist}{maximum pairwise distance between points (center-to-center, a.k.a. maximum spread, in original measurement units and MOA).}
\item{groupRect}{width and height of bounding box with diagonal and figure of merit FoM (average side length, in original measurement units and MOA).}
\item{groupRectMin}{width and height of minimum-area bounding box with diagonal and figure of merit FoM (average side length, in original measurement units and MOA).}
\item{minCircleRad}{radius for the minimum enclosing circle (in original measurement units and MOA).}
\item{confEll}{length of semi-major and semi-minor axis of the confidence ellipse (in original measurement units and MOA).}
\item{confEllRob}{length of semi-major and semi-minor axis of the confidence ellipse based on a robust estimate for the covariance matrix (in original measurement units and MOA).}
\item{confEllShape}{aspect ratio and flattening of the confidence ellipse.}
\item{confEllShapeRob}{aspect ratio and flattening of the confidence ellipse based on a robust estimate for the covariance matrix.}
\item{CEP}{estimate(s) for the circular error probable (CEP, in original measurement units and MOA).}

\item{ctr}{(x,y)-offset of group center relative to point of aim.}
\item{ctrXci}{parametric and bootstrap confidence intervals for center x-coordinate (1499 replicates).}
\item{ctrYci}{parametric and bootstrap confidence intervals for center y-coordinate (1499 replicates).}
\item{ctrRob}{robust estimate of group center offset relative to point of aim (MCD algorithm).}
\item{distPOA}{distance from group center to point of aim (in original measurement units and MOA).}
\item{distPOArob}{distance from robust estimate of group center to point of aim (in original measurement units and MOA).}
\item{Hotelling}{Hotelling's T^2-Test result from testing if group center equals point of aim.}
}
\seealso{
\code{\link[shotGroups]{groupShape}},
\code{\link[shotGroups]{groupLocation}},
\code{\link[shotGroups]{groupSpread}},
\code{\link[shotGroups]{compareGroups}},
\code{\link[shotGroups]{getDistToCtr}},
\code{\link[shotGroups]{getMaxPairDist}},
\code{\link[shotGroups]{getBoundingBox}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getMinCircle}},
\code{\link[shotGroups]{getConfEll}},
\code{\link[shotGroups]{getCEP}},
\code{\link[shotGroups]{getRayParam}},
\code{\link[shotGroups]{getMOA}},
\code{\link[graphics]{smoothScatter}},
\code{\link[mvoutlier]{chisq.plot}},
\code{\link[mvoutlier]{aq.plot}},
\code{\link[mvoutlier]{pcout}},
\code{\link[stats]{qqnorm}},
\code{\link[graphics]{hist}},
\code{\link[stats]{kernel}},
\code{\link[stats]{shapiro.test}},
\code{\link[energy]{mvnorm.etest}},
\code{\link[stats]{anova.mlm}},
\code{\link[boot]{boot}},
\code{\link[boot]{boot.ci}},
\code{\link[robustbase]{covMcd}}
}
\examples{
data(DFinch)

# select combined data from only first 2 series
DF  <- subset(DFinch, Series \%in\% 1:2)
res <- analyzeGroup(DF, conversion='yd2in')
names(res)
res$multNorm
res$corXY
res$ctrRob
res$ctrXci
res$ctrYci
}
