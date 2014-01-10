\name{DFsavage}
\alias{DFsavage}
\docType{data}
\title{Combined bullet hole data}
\description{
Example of a combined data frame from several original files, each with one group.
}
\usage{data(DFsavage)}
\format{
  A data frame with 180 observations on the following 10 variables.
  \describe{
    \item{\code{Point.X}}{a numerical vector of absolute x-coordinates of bullet holes.}
    \item{\code{Point.Y}}{a numerical vector of absolute y-coordinates of bullet holes.}
    \item{\code{Aim.X}}{a numerical vector of x-coordinates giving the point of aim.}
    \item{\code{Aim.Y}}{a numerical vector of y-coordinates giving the of point of aim.}
    \item{\code{Distance}}{a numerical vector giving the distance to the target.}
    \item{\code{Group}}{a factor with level \code{1}. This is the original \code{Group} variable as defined by OnTarget PC/TDS.}
    \item{\code{Bullet}}{a character vector describing the bullet type.}
    \item{\code{Origin}}{a factor with levels \code{1} \ldots \code{9}. This factor codes from which original output file the data is.}
    \item{\code{orgSer}}{a factor with levels \code{1.1} \ldots \code{9.1}. This factor results from \code{\link{droplevels}(\link{interaction}(Origin, Group))}, and codes each separate group across original files. The order of the factor levels is alphabetical.}
    \item{\code{Series}}{a factor with levels \code{1} \ldots \code{9}. This factor codes each separate group as defined by \code{orgSer}, but more conveniently as a number that runs consecutively across original files.}
 }
}
\details{
Several groups of shooting a Savage 12 FT/R rifle in .308 Win at distances from 100 to 300m. The measurement unit for coordinates is mm, for distance meters. Shots 1-5 in \code{Series} 4, and shots 1-3 in \code{Series} 7 moved the scope. \cr
This data frame is like those returned by \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}}, or \code{\link[shotGroups]{readDataMisc}} with option \code{combine=TRUE}. \cr
Data copyright Charles McMillan and Paul McMillan, 2008. \url{http://statshooting.com/}
}
\seealso{
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{analyzeGroup}},
\code{\link[shotGroups]{compareGroups}}
}
\examples{
data(DFsavage)
str(DFsavage)
}
\keyword{datasets}
