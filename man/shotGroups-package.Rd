\name{shotGroups-package}
\alias{shotGroups-package}
\alias{shotGroups}
\docType{package}
\title{Analyze shot group data}
\description{Analyze shot groups: functions to read in, plot, describe, and compare shot group data. This includes standard, but also nonparametric and robust methods statistical methods. Works with files exported by OnTarget PC and OnTarget TDS.}
\details{
\tabular{ll}{
Package: \tab shotGroups\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2012-02-21\cr
Depends: \tab R (>= 2.11), boot, KernSmooth, modeltools (>= 0.2-9)\cr
Imports: \tab coin, energy, graphics, grDevices, mvoutlier, robustbase, stats\cr
License: \tab GPL (>= 2)\cr
}
Use \code{help(package='shotGroups')} for a list of all functions. \cr \cr
If you do not yet have data imported into R:
\itemize{
\item if you have original output files from OnTarget PC/TDS: read in data from multiple files with \code{\link[shotGroups]{readDataOT1}} (for OnTarget PC v1.10), or with \code{\link[shotGroups]{readDataOT2}} (for OnTarget PC v2.10 and OnTarget TDS v3.71)
\item if you have other tab-delimited text files: read in data from multiple files with \code{\link[shotGroups]{readDataMisc}}
\item if your data is saved in some other text file format, consult \code{\link[utils]{read.table}}
\item combine data from different files with \code{\link[shotGroups]{combineData}}
}
If you already have full OnTarget PC/TDS data imported into an R:
\itemize{
\item make sure the data frame has these variables: \code{Group} (a \code{\link{factor}}), \code{Distance}, \code{Aim.X} and \code{Aim.Y} defining point of aim, \code{Point.X}, and \code{Point.Y} defining the bullet holes
\item perform a comprehensive numerical as well as graphical analysis of a group's shape, location, and spread with \code{\link[shotGroups]{analyzeGroup}}
\item numerically and visually compare different groups in terms of their shape, location, and spread with \code{\link[shotGroups]{compareGroups}}. This function requires an additional factor \code{Series} that identifies shot groups. OnTarget PC/TDS' \code{Group} variable identifies groups just within one file, \code{Series} is intended to number groups also across different original files. If you only have one original file, you can just copy variable \code{Groups} to \code{Series} in your data frame.
}
If you have data in the form of a matrix of (x,y)-coordinates relative to the point of aim, you can analyze different aspects of a group separately:
\itemize{
\item get numerical and graphical shape information on the distribution (find outliers and assess normality of distribution) with \code{\link[shotGroups]{groupShape}}
\item get location measures for the distribution with \code{\link[shotGroups]{groupLocation}}
\item get spread measures and their graphical representations with \code{\link[shotGroups]{groupSpread}}
}
}
\author{
Daniel Wollschlaeger

Maintainer: Daniel Wollschlaeger \email{dwoll@kuci.org}
}
\keyword{
package
}
\examples{
library(shotGroups)
example(groupShape)
example(groupLocation)
example(groupSpread)
example(analyzeGroup)
example(compareGroups)
}
