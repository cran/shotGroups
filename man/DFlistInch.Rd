\name{DFlistInch}
\alias{DFlistInch}
\docType{data}
\title{List containing several data frames with bullet hole data}
\description{
Example list containing several data frames with bullet hole data as produced by \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}}, or by \code{\link[shotGroups]{readDataMisc}}.
}
\usage{data(DFlistInch)}
\details{
Several groups of pistol shooting 9x19mm at 27yd The measurement unit for coordinates is inch, for distance yards. \cr \cr
This list can be used as an argument for \code{\link[shotGroups]{combineData}}.
}
\references{
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{readDataMisc}},
\code{\link[shotGroups]{readDataOT1}},
\code{\link[shotGroups]{readDataOT2}}
}
\examples{
data(DFlistInch)
str(DFlistInch)

## combine list of data frames to one single data frame
DFinch <- combineData(DFlistInch)
str(DFinch)
}
\keyword{datasets}
