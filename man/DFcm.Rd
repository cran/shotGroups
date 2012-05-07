\name{DFcm}
\alias{DFcm}
\docType{data}
\title{Combined bullet hole data}
\description{
Example of a combined data frame from several files exported by OnTarget PC/TDS, each with several groups of bullet holes.
}
\usage{data(DFcm)}
\format{
  A data frame with 487 observations on the following 13 variables.
  \describe{
    \item{\code{Project.Title} }{a character vector giving the OnTarget PC/TDS project title.}
    \item{\code{Group} }{a factor with levels \code{1} \code{2} \code{3}. This is the original \code{Group} variable as defined by OnTarget PC/TDS.}
    \item{\code{Ammunition} }{a character vector describing the ammo.}
    \item{\code{Distance} }{a numerical vector giving the distance to the target.}
    \item{\code{Aim.X} }{a numerical vector of x-coordinates giving the point of aim.}
    \item{\code{Aim.Y} }{a numerical vector of y-coordinates giving the of point of aim.}
    \item{\code{Center.X} }{a numerical vector of x-coordinates giving the group centers.}
    \item{\code{Center.Y} }{a numerical vector of y-coordinates giving the group centers.}
    \item{\code{Point.X} }{a numerical vector of absolute x-coordinates of bullet holes.}
    \item{\code{Point.Y} }{a numerical vector of absolute y-coordinates of bullet holes.}
    \item{\code{Origin} }{a factor with levels \code{1} \code{2} \code{3}. This factor codes from which original output file the data is.}
    \item{\code{orgSer} }{a factor with levels \code{1.1} \code{2.1} \code{3.1} \code{1.2} \code{2.2} \code{3.2} \code{1.3} \code{2.3} \code{3.3}. This factor results from \code{droplevels(interaction(Origin, Group))}, and codes each separate group across original files. The order of the factor levels is alphabetical.}
    \item{\code{Series} }{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9}. This factor codes each separate group as defined by \code{orgSer}, but more conveniently as a number that runs consecutively across original files.}
  }
}
\details{
Several groups of pistol shooting 9x19mm at 25m. The measurement unit for coordinates is cm, for distance meters. \cr
This data frame is the result from applying \code{\link[shotGroups]{combineData}} to a list of data frames as read in by \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}}, or by \code{\link[shotGroups]{readDataMisc}}.
}
\references{
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{analyzeGroup}},
\code{\link[shotGroups]{compareGroups}}
}
\examples{
data(DFcm)
str(DFcm)
}
\keyword{datasets}
