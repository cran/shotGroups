\name{DFtalon}
\alias{DFtalon}
\docType{data}
\title{Combined bullet hole data}
\description{
Example of a combined data frame from several original files, each with one group.
}
\usage{data(DFtalon)}
\format{
  A data frame with 180 observations on the following 10 variables.
  \describe{
    \item{\code{Point.X} }{a numerical vector of absolute x-coordinates of bullet holes.}
    \item{\code{Point.Y} }{a numerical vector of absolute y-coordinates of bullet holes.}
    \item{\code{Aim.X} }{a numerical vector of x-coordinates giving the point of aim.}
    \item{\code{Aim.Y} }{a numerical vector of y-coordinates giving the of point of aim.}
    \item{\code{Distance} }{a numerical vector giving the distance to the target.}
    \item{\code{Group} }{a factor with level \code{1}. This is the original \code{Group} variable as defined by OnTarget PC/TDS.}
    \item{\code{Bullet} }{a character vector describing the bullet type.}
    \item{\code{Origin} }{a factor with levels \code{1} \ldots \code{9}. This factor codes from which original output file the data is.}
    \item{\code{orgSer} }{a factor with levels \code{1.1} \ldots \code{9.1}. This factor results from \code{droplevels(interaction(Origin, Group))}, and codes each separate group across original files. The order of the factor levels is alphabetical.}
    \item{\code{Series} }{a factor with levels \code{1} \ldots \code{9}. This factor codes each separate group as defined by \code{orgSer}, but more conveniently as a number that runs consecutively across original files.}
  }
}
\details{
Several groups of air gun shooting at 10m. The measurement unit for coordinates is mm, for distance meters. \cr
This data frame is the result from applying \code{\link[shotGroups]{combineData}} to a list of data frames as read in by \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}}, or by \code{\link[shotGroups]{readDataMisc}}. \cr
Data copyright Charles McMillan and Paul McMillan, 2008. \url{http://statshooting.com/}
}
\references{
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{analyzeGroup}},
\code{\link[shotGroups]{compareGroups}}
}
\examples{
data(DFtalon)
str(DFtalon)
}
\keyword{datasets}
