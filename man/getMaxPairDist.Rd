\name{getMaxPairDist}
\alias{getMaxPairDist}
\title{Maximum pairwise distance for a set of points}
\description{Calculates the maximum of pairwise distances between points given a set of coordinates.}
\usage{getMaxPairDist(xy)}
\arguments{
  \item{xy }{a numerical (n x p)-matrix with the coordinates of n points in p-dimensional space (1 row of coordinates per point).}
}
\value{
A list with the following components
\item{d }{numerical value with the maximum pairwise distance between points.}
\item{idx }{a vector with two entries giving the row indices of the points that are farthest apart.}
}
\examples{
xy <- matrix(round(rnorm(20, 100, 15), 1), ncol=2)
(maxPD <- getMaxPairDist(xy))

plot(xy, asp=1, pch=16)
segments(x0=xy[maxPD$idx[1], 1], y0=xy[maxPD$idx[1], 2],
         x1=xy[maxPD$idx[2], 1], y1=xy[maxPD$idx[2], 2],
         col="green3", lwd=2)
}
