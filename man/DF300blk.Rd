\name{DF300BLK}
\alias{DF300BLK}
\docType{data}
\title{Combined bullet hole data}
\description{
Example of a data frame from one file with one group of bullet holes.}
\usage{data(DF300BLK)}
\format{
  A data frame with 20 observations on the following 9 variables.
  \describe{
    \item{\code{Group}}{a factor with level \code{1}.}
    \item{\code{Distance}}{a numerical vector giving the distance to the target.}
    \item{\code{Aim.X}}{a numerical vector of x-coordinates giving the point of aim.}
    \item{\code{Aim.Y}}{a numerical vector of y-coordinates giving the of point of aim.}
    \item{\code{Point.X}}{a numerical vector of absolute x-coordinates of bullet holes.}
    \item{\code{Point.Y}}{a numerical vector of absolute y-coordinates of bullet holes.}
    \item{\code{Origin}}{a factor with level \code{1}. This factor codes from which original output file the data is.}
    \item{\code{orgSer}}{a factor with level \code{1.1}. This factor results from \code{\link{droplevels}(\link{interaction}(Origin, Group))}, and codes each separate group across original files. The order of the factor levels is alphabetical.}
    \item{\code{Series}}{a factor with level \code{1}. This factor codes each separate group as defined by \code{orgSer}, but more conveniently as a number that runs consecutively across original files.}
  }
}
\details{
One group of shooting a Noveske AR-15 rifle in 300BLK at 100yd with factory ammunition. The measurement unit for coordinates is inch, for distance yards. \cr
This data frame is like those returned by \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}}, or \code{\link[shotGroups]{readDataMisc}} with option \code{combine=TRUE}. \cr
Data courtesy of David Bookstaber, 2013. \url{http://ballistipedia.com/}
}
\seealso{
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{analyzeGroup}},
\code{\link[shotGroups]{compareGroups}}
}
\examples{
data(DF300BLK)
str(DF300BLK)
}
\keyword{datasets}
