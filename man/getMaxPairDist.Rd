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
getMaxPairDist(xy)
}
