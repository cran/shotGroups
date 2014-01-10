\name{combineData}
\alias{combineData}
\title{Combine list of data frames into one}
\description{Combines a list of data frames (the result from using \code{\link[shotGroups]{readDataOT1}}, \code{\link[shotGroups]{readDataOT2}}, or \code{\link[shotGroups]{readDataMisc}} with \code{combine=FALSE}) into one big data frame.}
\usage{combineData(DFs)}
\arguments{
  \item{DFs}{ a list of data frames with a shared set of variables.}
}
\details{Assumes that the data frames in the list have a non-empty set of shared variables. Among them at least either \code{Point.X}, \code{Point.Y} or \code{X}, \code{Y} defining the bullet holes. To be useful for functions \code{\link[shotGroups]{analyzeGroup}} or \code{\link[shotGroups]{compareGroups}}, the data frames also need variables \code{Group}, \code{Distance}, \code{Aim.X}, \code{Aim.Y} defining point of aim. If they are missing, they are set to defaults (1, 100, (0,0), respectively).}
\value{
A data frame with the shared set of variables. In addition, it also contains factors identifying the original file (\code{Origin}), and a single consecutive numbering of all groups from different data frames (\code{Series}).
\item{\ldots}{the shared set of variables from the the data frames in the list.}
\item{Group}{a factor that is the original \code{Group} variable as defined by OnTarget PC/TDS.}
\item{Origin}{a factor that codes from which original file the data is.}
\item{orgSer}{a factor that results from \code{\link{droplevels}(\link{interaction}(Origin, Group))}, and codes each separate group across original files. The order of the factor levels is alphabetical.}
\item{Series}{a factor that codes each separate group as defined by \code{orgSer}, but more conveniently as a number that runs consecutively across original files.}
}
\seealso{
\code{\link[shotGroups]{readDataMisc}},
\code{\link[shotGroups]{readDataOT1}},
\code{\link[shotGroups]{readDataOT2}},
\code{\link[shotGroups]{analyzeGroup}},
\code{\link[shotGroups]{compareGroups}}
}
\examples{
## combine list of data frames to one single data frame
data(DFlistCm)
DFcm <- combineData(DFlistCm)
str(DFcm)
head(DFcm)
}
