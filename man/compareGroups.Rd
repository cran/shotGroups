\name{compareGroups}
\alias{compareGroups}
\title{Compare bullet hole groups}
\description{Numerically and graphically compares the distributions of up to 15 groups of bullet holes.}
\usage{compareGroups(DF, plots = TRUE, xyTopLeft = TRUE,
              ABalt = c('two.sided', 'less', 'greater'),
              Walt = c('two.sided', 'less', 'greater'),
              conversion = 'm2cm')}
\arguments{
  \item{DF }{a data frame containing (at least) these variables: \code{Series} (a \code{\link{factor}}), \code{Distance}, \code{Aim.X}, \code{Aim.Y}defining point of aim, \code{Point.X}, and \code{Point.Y} defining the bullet holes.}
  \item{plots }{a logical value indicating whether diagrams should be shown.}
  \item{xyTopLeft }{a logical value indicating whether the origin of the absolute coordinate system is in the top-left corner. This is the default for data exported by OnTarget PC/TDS.}
  \item{ABalt }{a character string indicating the hypothesis for the Ansari-Bradley-Test for equal variances. Only used when exactly 2 groups are compared.}
  \item{Walt }{a character string indicating the hypothesis for the Wilcoxon-Rank-Sum-Test for equality of average distance to group center (equivalent to the Mann-Whitney-U-Test). Only used when exactly 2 groups are compared.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
}
\details{
OnTarget PC/TDS' \code{Group} variable identifies groups just within one file, whereas factor \code{Series} is taken to number groups also across different original files. If you only have one original file, you can just copy variable \code{Groups} to \code{Series} in your data frame. \cr \cr
In addition to the numerical results listed below, this function produces the following diagrams:
\itemize{
\item a scatterplot showing all groups, as well as their respective center and characteristic ellipse
\item a scatterplot showing all groups, as well as their respective minimum bounding box, minimum enclosing circle, and maximum group spread
\item a scatterplot showing all groups, as well as their respective minimum enclosing circle and circle with average distance to center
}
}
\value{
A list with the results from numerical comparisons and statistical tests.
\item{ctr }{group centers relative to the respective point of aim.}
\item{distPOA }{distances from group centers to point of aim in original measurement units and MOA.}
\item{MANOVA }{MANOVA result from testing equality of group centers (test statistic is Wilk's lambda).}
\item{sdXY }{group standard deviations for the x- and y-coordinates.}
\item{covXY }{group covariance matrices for the (x,y)-coordinates.}
\item{meanDistToCtr }{average distances from points to their group centers in original measurement units and MOA.}
\item{maxPairDist }{maximum pairwise distance between points for each group in original measurement units and MOA.}
\item{groupWidth }{width of minimum bounding box for each group in original measurement units and MOA.}
\item{groupHeight }{height of minimum bounding box for each group in original measurement units and MOA.}
\item{minCircleRad }{radius for the minimum enclosing circle for each group in original measurement units and MOA.}
\item{AnsariX }{Ansari-Bradley-Test result from testing equality of group variances for x-coordinates. When two groups are compared.}
\item{AnsariY }{Ansari-Bradley-Test result from testing equality of group variances for y-coordinates. When two groups are compared.}
\item{Wilcoxon }{Wilcoxon-Rank-Sum-Test result from testing equality of average point distances to their respective group center. When two groups are compared.}
\item{FlignerX }{Fligner-Killeen-Test result from testing equality of group variances for x-coordinates. When more than two groups are compared.}
\item{FlignerY }{Fligner-Killeen-Test result from testing equality of group variances for y-coordinates. When more than two groups are compared.}
\item{Kruskal }{Kruskal-Wallis-Test result from testing equality of average point distances to their respective group center. When more than two groups are compared.}
}
\seealso{
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{analyzeGroup}},
\code{\link[shotGroups]{getDistToCtr}},
\code{\link[shotGroups]{getMaxPairDist}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getMinCircle}},
\code{\link[shotGroups]{getMOA}},
\code{\link[shotGroups]{drawEllipse}},
\code{\link[stats]{anova.mlm}},
\code{\link[coin]{ansari_test}},
\code{\link[coin]{fligner_test}},
\code{\link[coin]{wilcox_test}},
\code{\link[coin]{kruskal_test}}
}
\examples{
data(DFinch)

# select only first 2 series
DFsub <- subset(DFinch, Series \%in\% c(2:3))
cmp   <- compareGroups(DFsub, conversion='yd2in')
names(cmp)
cmp$ctr
cmp$meanDistToCtr
cmp$Wilcoxon
}
