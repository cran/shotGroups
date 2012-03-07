\name{getDistToCtr}
\alias{getDistToCtr}
\title{Distances to center for a set of points}
\description{Calculates the respective distances for a set of points to their center.}
\usage{getDistToCtr(xy)}
\arguments{
  \item{xy }{a numerical (n x p)-matrix with the coordinates of n points in p-dimensional space (1 row of coordinates per point).}
}
\value{A numerical vector with the distances from each point to the center of the set.}
\examples{
xy <- matrix(round(rnorm(20, 100, 15), 1), ncol=2)
getDistToCtr(xy)
}
