\name{DFlistTalon}
\alias{DFlistTalon}
\docType{data}
\title{List containing several data frames with bullet hole data}
\description{
Example list containing several data frames with bullet hole data as produced by \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}}, or by \code{\link[shotGroups]{readDataMisc}}.
}
\usage{data(DFlistTalon)}
\details{
Several groups of air gun shooting at 10m. The measurement unit for coordinates is mm, for distance meters. \cr
This list can be used as an argument for \code{\link[shotGroups]{combineData}}. \cr
Data copyright Charles McMillan and Paul McMillan, 2008. \url{http://statshooting.com/}
}
\references{
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{readDataMisc}},
\code{\link[shotGroups]{readDataOT1}},
\code{\link[shotGroups]{readDataOT2}}
}
\examples{
data(DFlistTalon)
str(DFlistTalon)

## combine list of data frames to one single data frame
DFtalon <- combineData(DFlistTalon)
str(DFtalon)
}
\keyword{datasets}
