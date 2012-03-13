\name{analyzeGroup}
\alias{analyzeGroup}
\title{Analysis for a single group of bullet holes}
\description{Performs a comprehensive numerical and graphical analysis of a single group of bullet holes.}
\usage{analyzeGroup(DF, xyTopLeft = TRUE, conversion = 'm2cm', bandW = 0.5)}
\arguments{
  \item{DF }{a data frame containing (at least) these variables: \code{Distance}, \code{Aim.X}, \code{Aim.Y} defining point of aim, \code{Point.X}, and \code{Point.Y} defining the bullet holes.}
  \item{xyTopLeft }{a logical value indicating whether the origin of the absolute coordinate system is in the top-left corner. This is the default for data exported by OnTarget PC/TDS.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{bandW }{for argument \code{bandwith} of \code{\link[graphics]{smoothScatter}}.}
}
\details{
Robust estimates for the group center and the covariance matrix of (x,y)-coordinates are from \code{\link[robustbase]{covMcd}} using the MCD algorithm. \cr
This function is a wrapper for \code{\link[shotGroups]{groupShape}}, \code{\link[shotGroups]{groupLocation}}, and \code{\link[shotGroups]{groupSpread}}. \cr
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
\item a histogram for distances to group center including a nonparametric kernel density estimate
}
}
\value{
A list with the results from the numerical analyses and statistical tests.
\item{corXY }{correlation matrix of (x,y)-coordinates.}
\item{corXYrob }{robust estimate of correlation matrix of (x,y)-coordinates.}
\item{Outliers }{a vector of row indices for observations identified as outliers.}
\item{ShapiroX }{Shapiro-Wilk-Test result for normality of x-coordinates.}
\item{ShapiroY }{Shapiro-Wilk-Test result for normality of y-coordinates.}
\item{multNorm }{E-statistic-Test result for multivariate normality of (x,y)-coordinates.}
\item{sdXY }{standard deviations of x- and y coordinates (in original measurement units and MOA).}
\item{CIsdXY }{parametric 95\%-confidence interval for the standard deviations of x- and y coordinates (in original measurement units).}
\item{sdXYrob }{robust standard deviations of x- and y-coordinates (in original measurement units and MOA).}
\item{covXY }{covariance matrix of (x,y)-coordinates.}
\item{covXYrob }{robust estimate of covariance matrix of (x,y)-coordinates.}
\item{meanDistToCtr }{average distance from points to their center (in original measurement units and MOA).}
\item{maxPairDist }{maximum pairwise distance between points (center-to-center, a.k.a. maximum spread, in original measurement units and MOA).}
\item{groupRect }{width and height of minimum bounding box (in original measurement units and MOA).}
\item{minCircleRad }{radius for the minimum enclosing circle (in original measurement units and MOA).}
\item{confEll }{length of semi-major and semi-minor axis of the 50\%-confidence ellipse (in original measurement units and MOA).}
\item{confEllRob }{length of semi-major and semi-minor axis of the 50\%-confidence ellipse based on a robust estimate for the covariance matrix (in original measurement units and MOA).}
\item{confEllShape }{aspect ratio and flattening of the confidence ellipse.}
\item{confEllShapeRob }{aspect ratio and flattening of the confidence ellipse based on a robust estimate for the covariance matrix.}
\item{CEPrand }{modified RAND R-234 estimate for the circular error probable (CEP, in original measurement units and MOA).}
\item{ctr }{(x,y)-offset of group center relative to point of aim.}
\item{ctrRob }{robust estimate of group center offset relative to point of aim.}
\item{distPOA }{distance from group center to point of aim (in original measurement units and MOA).}
\item{distPOArob }{distance from robust estimate of group center to point of aim (in original measurement units and MOA).}
\item{Hotelling }{Hotelling's T^2-Test result from testing if group center equals point of aim.}
\item{CItX }{95\% t-confidence interval for x-coordinates.}
\item{CItY }{95\% t-confidence interval for y-coordinates.}
\item{CIbootX }{95\% bootstrap-confidence intervals for x-coordinates (percentile and BCa, 1499 replicates).}
\item{CIbootY }{95\% bootstrap-confidence intervals for y-coordinates (percentile and BCa, 1499 replicates).}
}
\seealso{
\code{\link[shotGroups]{groupShape}},
\code{\link[shotGroups]{groupLocation}},
\code{\link[shotGroups]{groupSpread}},
\code{\link[shotGroups]{compareGroups}},
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{getDistToCtr}},
\code{\link[shotGroups]{getMaxPairDist}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getMinCircle}},
\code{\link[shotGroups]{getConfEll}},
\code{\link[shotGroups]{getCEP}},
\code{\link[shotGroups]{getMOA}},
\code{\link[graphics]{smoothScatter}},
\code{\link[mvoutlier]{aq.plot}},
\code{\link[mvoutlier]{chisq.plot}},
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

# select combined data from only first 4 series
DF  <- subset(DFinch, Series \%in\% c(1:4))
res <- analyzeGroup(DF, conversion='yd2in')
names(res)
res$multNorm
res$corXY
res$ctrRob
res$CIbootX
res$CIbootY
}
