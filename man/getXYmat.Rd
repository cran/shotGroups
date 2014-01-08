\name{getXYmat}
\alias{getXYmat}
\title{Extract (x,y)-coordinates (relative to point of aim) from a data frame}
\description{Extracts (x,y)-coordinates of the points of impact (relative to the point of aim) from a data frame and returns them as a matrix.}
\usage{getXYmat(DF, xyTopLeft = TRUE, relPOA = TRUE)}
\arguments{
  \item{DF}{a data frame containing (at least) either the variables \code{Point.X}, \code{Point.Y} or \code{X}, \code{Y} defining the bullet holes. Variables \code{Aim.X}, \code{Aim.Y} can define the point of aim. If missing, point of aim is assumed to be in (0,0).}
  \item{xyTopLeft}{a logical value indicating whether the origin of the absolute coordinate system is in the top-left corner. This is the default for data exported by OnTarget PC/TDS.}
  \item{relPOA}{a logical value indicating whether returned coordinates should be relative to the point of aim.}
}
\seealso{
\code{\link[shotGroups]{groupLocation}},
\code{\link[shotGroups]{groupShape}},
\code{\link[shotGroups]{groupSpread}}
}
\value{A numerical matrix with the (x,y)-coordinates.}
\examples{
data(DFcm)

# select data from only first series
DFsub <- subset(DFcm, Series == 1)
getXYmat(DFsub, xyTopLeft=TRUE, relPOA=TRUE)
}
