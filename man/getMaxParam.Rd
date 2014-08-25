\name{getMaxParam}
\alias{getMaxParam}
\alias{getMaxParam.data.frame}
\alias{getMaxParam.default}
\title{Estimate Maxwell-Boltzmann parameters sigma, mean and standard deviation}
\description{Estimates the scale parameter sigma of the Maxwell-Boltzmann distribution together with the radial mean MR and radial standard deviation RSD, including parametric confidence intervals.}
\usage{getMaxParam(xyz, level = 0.95, mu, doRob = FALSE)

\method{getMaxParam}{data.frame}(xyz, level = 0.95, mu, doRob = FALSE)

\method{getMaxParam}{default}(xyz, level = 0.95, mu, doRob = FALSE)}
\arguments{
  \item{xyz}{either a numerical (n x 3)-matrix with the coordinates of n points (1 row of (x,y,z)-coordinates per point), or a data frame with either the variables \code{X}, \code{Y}, \code{Z} or \code{Point.X}, \code{Point.Y}, \code{Point.Z}.}
  \item{level}{a numerical value with the coverage for the confidence intervals for sigma, MR, RSD.}
  \item{mu}{true distribution center (optional). See details.}
  \item{doRob}{logical: use robust estimation of covariance matrix as basis for estimators?}
}
\details{When the true mean \code{mu} of the distribution is given, the sigma estimate uses the sum of squared radii for the variance estimate (the total un-corrected variance of the coordinates), and employs the c4(3*N+1) correction factor and for taking the square root. When \code{mu} is missing, the sum of squared radii is Bessel-corrected for estimating the center, and c4(3*N-2) is used.

MR and RSD each are a constant fraction of sigma.

The robust estimate for the covariance matrix of (x,y,z)-coordinates is from \code{\link[robustbase]{covMcd}} using the MCD algorithm.
}
\value{
A list with the estimates for sigma, RSD, and MR including the confidence intervals.
\item{sigma}{A vector with the sigma estimate and confidence interval bounds as named elements \code{sigma}, \code{sigCIlo}, \code{sigCIup}.}
\item{MR}{A vector with the MR estimate and confidence interval bounds as named elements \code{MR}, \code{MRciLo}, \code{MRciUp}.}
\item{RSD}{A vector with the RSD estimate and confidence interval bounds as named elements \code{RSD}, \code{RSDciLo}, \code{RSDciUp}.}
}
\references{
\url{http://reference.wolfram.com/language/ref/MaxwellDistribution.html}
}
\seealso{
\code{\link[shotGroups]{Maxwell}},
\code{\link[shotGroups]{getCEP}},
\code{\link[shotGroups]{getHitProb}},
\code{\link[robustbase]{covMcd}}
}
\examples{
# coordinates given by a matrix
xyz <- matrix(round(rnorm(60, 0, 5), 2), ncol=3)
getMaxParam(xyz, doRob=FALSE)
}
