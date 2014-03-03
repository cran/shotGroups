\name{getRayParam}
\alias{getRayParam}
\alias{getRayParam.data.frame}
\alias{getRayParam.default}
\title{Estimate Rayleigh parameters sigma, RSD, and MR}
\description{Estimates the radial precision parameter sigma of the Rayleigh distribution together with the Rayleigh radial standard deviation RSD and the Rayleigh mean MR, including parametric confidence intervals.}
\usage{getRayParam(xy, level = 0.95, accuracy = FALSE)

\method{getRayParam}{data.frame}(xy, level = 0.95, accuracy = FALSE)

\method{getRayParam}{default}(xy, level = 0.95, accuracy = FALSE)}
\arguments{
  \item{xy}{either a numerical (n x p)-matrix with the coordinates of n points (1 row of p coordinates per point), or a data frame with either the variables \code{X}, \code{Y} or \code{Point.X}, \code{Point.Y} (optionally \code{Z} or \code{Point.Z}).}
  \item{level}{a numerical value with the coverage for the confidence intervals for sigma, RSD, MR.}
  \item{accuracy}{logical: take systematic location bias into account for the sigma estimate? See details.}
}
\details{By default, the reported sigma estimate does not take into account accuracy, i.e., any systematic location bias. The true mean of the distribution is then assumed to coincide with the group center. With \code{accuracy=TRUE}, the true center is set to be the point of aim in the origin 0. RSD and MR each are a constant fraction of sigma.

The sigma estimate uses the sum of squared radii as a variance estimate (equal to the total un-corrected variance of the coordinates), and then employs a correction factor for taking the square root. With \code{accuracy=FALSE}, MR is the distance to the group center, otherwise MR is the distance to the point of aim in the origin 0.
}
\value{
A list with the estimates for sigma, RSD, and MR including the confidence intervals.
\item{sigma}{A vector with the sigma estimate and confidence interval bounds as named elements \code{sigma}, \code{sigCIlo}, \code{sigCIup}.}
\item{RSD}{A vector with the RSD estimate and confidence interval bounds as named elements \code{RSD}, \code{RSDciLo}, \code{RSDciUp}.}
\item{MR}{A vector with the MR estimate and confidence interval bounds as named elements \code{MR}, \code{MRciLo}, \code{MRciUp}.}
}
\references{
\url{http://ballistipedia.com/index.php?title=Measuring_Precision}
}
\seealso{
\code{\link[shotGroups]{getCEP}},
\code{\link[shotGroups]{getHitProb}},
\code{\link[shotGroups]{groupSpread}}
}
\examples{
# coordinates given by a suitable data frame
getRayParam(DFtalon, level=0.95, accuracy=FALSE)

# coordinates given by a matrix
\dontrun{
xy <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
getRayParam(xy, level=0.95, accuracy=FALSE)
}
}
