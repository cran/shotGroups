\name{DFlistSavage}
\alias{DFlistSavage}
\docType{data}
\title{List containing several data frames with bullet hole data}
\description{
Example list containing several data frames with bullet hole data as produced by \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}}, or by \code{\link[shotGroups]{readDataMisc}}.
}
\usage{data(DFlistSavage)}
\details{
Several groups of rifle shooting at distances from 100m to 300m. The measurement unit for coordinates is mm, for distance meters. Shots 1-5 in \code{file4}, and shots 1-3 in \code{file7} moved the scope. \cr
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
data(DFlistSavage)
str(DFlistSavage)

## combine list of data frames to one single data frame
DFsavage <- combineData(DFlistSavage)
str(DFsavage)
}
\keyword{datasets}
