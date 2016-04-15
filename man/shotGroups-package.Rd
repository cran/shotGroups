\name{shotGroups-package}
\alias{shotGroups-package}
\alias{shotGroups}
\docType{package}
\title{Analyze shot group data}
\description{The \code{shotGroups} package provides functions to read in, plot, statistically describe, analyze, and compare shooting data with respect to group shape, precision, and accuracy. This includes graphical methods, descriptive statistics, and inference tests using standard, but also non-parametric and robust statistical methods. The data can be imported from files produced by OnTarget PC and OnTarget TDS or Taran, or from custom data files in text format with a similar structure. Supports inference from range statistics. Includes web-based graphical user interface.}
\details{
\tabular{ll}{
Package: \tab shotGroups\cr
Type: \tab Package\cr
Version: \tab 0.7.1\cr
Date: \tab 2016-04-25\cr
Depends: \tab R (>= 3.1.1)\cr
Imports: \tab boot, coin, CompQuadForm, graphics, grDevices, KernSmooth, robustbase, stats, tools, utils\cr
Suggests: \tab knitr, energy, mvoutlier, shiny\cr
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
groupSpread(DFcciHV, dstTarget=100, conversion='yd2in', bootCI='none')
}
