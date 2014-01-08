\name{compareGroups}
\alias{compareGroups}
\title{Compare bullet hole groups}
\description{Numerically and graphically compare accuracy, precision, and distribution shape of up to 15 groups of bullet holes.}
\usage{compareGroups(DF, plots = TRUE, xyTopLeft = TRUE,
              ABalt = c('two.sided', 'less', 'greater'),
              Walt = c('two.sided', 'less', 'greater'),
              conversion = 'm2cm')}
\arguments{
  \item{DF}{a data frame containing (at least) these variables: \code{Series} (a \code{\link{factor}}), and either \code{Point.X}, \code{Point.Y} or \code{X}, \code{Y} defining the bullet holes. Variables \code{Distance} (distance to target), \code{Aim.X}, \code{Aim.Y} (point of aim) are useful - if they are missing, a warning is given and a default assumed.}
  \item{plots}{a logical value indicating whether diagrams should be shown.}
  \item{xyTopLeft}{a logical value indicating whether the origin of the absolute coordinate system is in the top-left corner. This is the default for data exported by OnTarget PC/TDS.}
  \item{ABalt}{a character string indicating the hypothesis for the Ansari-Bradley-Test for equal variances. Only used when exactly 2 groups are compared.}
  \item{Walt}{a character string indicating the hypothesis for the Wilcoxon-Rank-Sum-Test for equality of average distance to group center (equivalent to the Mann-Whitney-U-Test). Only used when exactly 2 groups are compared.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
}
\details{
OnTarget PC/TDS' \code{Group} variable identifies groups just within one file, whereas factor \code{Series} is taken to number groups also across different original files. If your data was read with \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}} or \code{\link[shotGroups]{readDataMisc}}, \code{Series} is added automatically. Otherwise, you can just copy variable \code{Groups} to \code{Series} in a data frame called \code{shots} with \code{shots$Series} \code{<-} \code{shots$Group}.\cr
If the data is missing information about the point of aim, (0,0) is assumed. If distance to target is missing, 100 is assumed.\cr \cr
In addition to the numerical results listed below, this function produces the following diagrams:
\itemize{
\item a scatterplot showing all groups as well as their respective center and 50\%-confidence ellipse
\item a scatterplot showing all groups as well as their respective minimum bounding box, minimum enclosing circle, and maximum group spread
\item a scatterplot showing all groups as well as their respective minimum enclosing circle and circle with average distance to center
}
}
\value{
A list with the results from numerical comparisons and statistical tests.
\item{ctr}{group centers relative to their respective point of aim.}
\item{distPOA}{distances from group centers to point of aim (in original measurement units and MOA).}
\item{MANOVA}{MANOVA result from testing equality of group centers relative to their respective point of aim (test statistic is Wilk's lambda).}
\item{corXY}{group correlation matrices for the (x,y)-coordinates.}
\item{sdXY}{group standard deviations for the x- and y-coordinates (in original measurement units).}
\item{sdXYmoa}{group standard deviations for the x- and y-coordinates (in MOA).}
\item{sdXYci}{group parametric (chi^2) 95\%-confidence intervals for the standard deviations of x- and y coordinates (in original measurement units).}
\item{sdXYciMOA}{group parametric (chi^2) 95\%-confidence intervals for the standard deviations of x- and y coordinates (in MOA).}
\item{meanDistToCtr}{average distances from points to their respective group center (in original measurement units and MOA).}
\item{maxPairDist}{maximum pairwise distance between points for each group (center-to-center, a.k.a. maximum spread, in original measurement units and MOA).}
\item{bbFoM}{minimum-area bounding box figure of merit (average side length) for each group (in original measurement units and MOA).}
\item{bbDiag}{minimum-area bounding box diagonal length for each group (in original measurement units and MOA).}
\item{minCircleRad}{radius of the minimum enclosing circle for each group (in original measurement units and MOA).}
\item{CEPray}{Rayleigh estimate for the 50\% circular error probable (CEP) in each group (in original measurement units and MOA).}
\item{AnsariX}{Ansari-Bradley-Test result from testing equality of group variances for x-coordinates. When two groups are compared.}
\item{AnsariY}{Ansari-Bradley-Test result from testing equality of group variances for y-coordinates. When two groups are compared.}
\item{Wilcoxon}{Wilcoxon-Rank-Sum-Test result from testing equality of average point distances to their respective group center. When two groups are compared.}
\item{FlignerX}{Fligner-Killeen-Test result from testing equality of group variances for x-coordinates. When more than two groups are compared.}
\item{FlignerY}{Fligner-Killeen-Test result from testing equality of group variances for y-coordinates. When more than two groups are compared.}
\item{Kruskal}{Kruskal-Wallis-Test result from testing equality of average point distances to their respective group center. When more than two groups are compared.}
}
\seealso{
\code{\link[shotGroups]{analyzeGroup}},
\code{\link[shotGroups]{getDistToCtr}},
\code{\link[shotGroups]{getMaxPairDist}},
\code{\link[shotGroups]{getMinBBox}},
\code{\link[shotGroups]{getMinCircle}},
\code{\link[shotGroups]{getCEP}},
\code{\link[shotGroups]{getMOA}},
\code{\link[shotGroups]{getRayParam}},
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
DFsub <- subset(DFinch, Series \%in\% 2:3)
cmp   <- compareGroups(DFsub, conversion='yd2in')
names(cmp)
cmp$ctr
cmp$meanDistToCtr
cmp$CEPray
cmp$Wilcoxon
}
