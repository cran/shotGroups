\name{getXYmat}
\alias{getXYmat}
\title{Matrix of (x,y)-coordinates relative to point of aim}
\description{Returns the matrix of (x,y)-coordinates of the points of impact relative to the point of aim. This matrix can then be used in functions like \code{\link[shotGroups]{groupLocation}}, \code{\link[shotGroups]{groupShape}}, or \code{\link[shotGroups]{groupSpread}}.}
\usage{getXYmat(DF, xyTopLeft = TRUE)}
\arguments{
  \item{DF }{a data frame containing (at least) these variables: \code{Aim.X}, \code{Aim.Y} defining point of aim, \code{Point.X}, and \code{Point.Y} defining the bullet holes.}
  \item{xyTopLeft }{a logical value indicating whether the origin of the absolute coordinate system is in the top-left corner. This is the default for data exported by OnTarget PC/TDS.}
}
\seealso{
\code{\link[shotGroups]{groupLocation}},
\code{\link[shotGroups]{groupShape}},
\code{\link[shotGroups]{groupSpread}}
}
\value{A numerical matrix with the (x,y)-coordinates.}
\examples{
data(DFcm)

# select combined data from only first 2 series
DF <- subset(DFcm, Series \%in\% 1:2)
xy <- getXYmat(DF, xyTopLeft=TRUE)
groupLocation(xy, plots=2)
}
