\name{getRayParam}
\alias{getRayParam}
\alias{getRayParam.data.frame}
\alias{getRayParam.default}
\title{Estimate Rayleigh parameters sigma, RSD, and MR}
\description{Estimates the radial precision parameter sigma of the Rayleigh distribution together with the Rayleigh radial standard deviation RSD and the Rayleigh mean MR, including parametric confidence intervals.}
\usage{getRayParam(xy, type = c('Rayleigh', 'Gauss'), level = 0.95, accuracy = FALSE)

\method{getRayParam}{data.frame}(xy, type = c('Rayleigh', 'Gauss'), level = 0.95, accuracy = FALSE)

\method{getRayParam}{default}(xy, type = c('Rayleigh', 'Gauss'), level = 0.95, accuracy = FALSE)}
\arguments{
  \item{xy}{either a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point), or a data frame with either the variables \code{X}, \code{Y} or \code{Point.X}, \code{Point.Y}.}
  \item{type}{string vector indicating which sigma estimate to report (see details).}
  \item{level}{a numerical value with the coverage for the confidence intervals for sigma, RSD, MR.}
  \item{accuracy}{a logical value indicating whether systematic location bias should be taken into account with the Rayleigh sigma estimate (see details).}
}
\details{By default, the reported sigma estimate does not take into account accuracy, i.e., any systematic location bias. The Rayleigh estimate combined with \code{accuracy=TRUE} can incorporate systematic deviation from point of aim. RSD and MR each are a constant fraction of sigma.
\itemize{
\item The Rayleigh estimate uses the sum of squared radii as a variance estimate, and then employs a correction factor for taking the square root. With \code{accuracy=FALSE}, the data is centered first, otherwise the radius is the distance to (0,0), assumed to be the point of aim.
\item The Gauss estimate uses the average x- and y-variance as a variance estimate, and then employs a correction factor for taking the square root. The option \code{accuracy} is ignored, the behavior is always as \code{accuracy=FALSE}.
}
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
\code{\link[shotGroups]{groupSpread}}
}
\examples{
# coordinates given by a suitable data frame
getRayParam(DFtalon, type='Gauss', level=0.95, accuracy=FALSE)

# coordinates given by a matrix
\dontrun{
xy <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
getRayParam(xy, type='Gauss', level=0.95, accuracy=FALSE)
}
}
