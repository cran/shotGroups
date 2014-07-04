\name{shotGroups-package}
\alias{shotGroups-package}
\alias{shotGroups}
\docType{package}
\title{Analyze shot group data}
\description{The \code{shotGroups} package provides functions to read in, plot, statistically describe, analyze, and compare shooting data with respect to group shape, precision, and accuracy. This includes graphical methods, descriptive statistics, and inference tests using standard, but also nonparametric and robust statistical methods. The data can be imported from files produced by OnTarget PC and OnTarget TDS or from custom data files in text format with a similar structure.}
\details{
\tabular{ll}{
Package: \tab shotGroups\cr
Type: \tab Package\cr
Version: \tab 0.5\cr
Date: \tab 2014-03-03\cr
Depends: \tab R (>= 2.15.2)\cr
Imports: \tab boot, coin, CompQuadForm, energy, graphics, grDevices, mvoutlier, robustbase, stats, KernSmooth\cr
License: \tab GPL (>= 2)\cr
}
Use \code{help(package='shotGroups')} for a list of all functions and links to the detailed help pages with information on options, usage and output. For further explanantions and an example walkthrough, see \code{vignette('shotGroups')}.
}
\author{
Daniel Wollschlaeger

Maintainer: Daniel Wollschlaeger \email{dwoll@kuci.org}
}
\keyword{
package
}
\examples{
groupSpread(DFcciHV, dstTarget=100, conversion='yd2in', bootCI='basic')
}
