\name{getHitProb}
\alias{getHitProb}
\alias{getHitProb.data.frame}
\alias{getHitProb.default}
\title{Hit probability within given region}
\description{Calculates the hit probability within a circular or spherical region for a given group of two- or three-dimensional coordinates.}
\usage{
getHitProb(xy, r=1, unit = 'unit', dstTarget = 100, conversion = 'm2cm',
           accuracy = FALSE, type = 'CorrNormal', doRob = FALSE)

\method{getHitProb}{data.frame}(xy, r=1, unit = 'unit', dstTarget = 100, conversion = 'm2cm',
           accuracy = FALSE, type = 'CorrNormal', doRob = FALSE)

\method{getHitProb}{default}(xy, r=1, unit = 'unit', dstTarget = 100, conversion = 'm2cm',
           accuracy = FALSE, type = 'CorrNormal', doRob = FALSE)
}
\arguments{
  \item{xy}{either a numerical (n x p)-matrix with the coordinates of n points (1 row of p coordinates per point), or a data frame with either the variables \code{X}, \code{Y} or \code{Point.X}, \code{Point.Y} (optionally \code{Z} or \code{Point.Z}).}
  \item{r}{a numerical value with the radius for the circle/sphere that defines the region for which the probability should be calculated.}
  \item{unit}{Measurement unit for radius \code{r}. Default \code{'unit'} indicates that the measurement unit for (x,y)-coordinates given in \code{conversion}. Possible values are \code{'unit', 'm', 'cm', 'mm', 'yd', 'ft', 'in', 'MOA', 'SMOA', 'milrad'}.}
  \item{dstTarget}{a numerical value with the distance to the target - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the shot-coordinates - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{accuracy}{logical: take systematic location bias into account? See details.}
  \item{type}{string vector indicating which estimation method to use. Possible values are \code{'CorrNormal', 'GrubbsPearson', 'GrubbsPatnaik', 'GrubbsLiu', 'Rayleigh'}.}
  \item{doRob}{logical: use robust estimation of center and covariance matrix as basis for estimates?}
}
\details{
The estimated hit probability is the same as the estimated fraction of shots with a distance to their group center of at most \code{r}.

By default, the estimated hit probability does not take into account accuracy, i.e., any systematic location bias. It is then based on the centered data. Set \code{accuracy=TRUE} to incorporate systematic deviation from point of aim.
\itemize{
\item \code{CorrNormal}: For \code{accuracy=FALSE} and two-dimensional data, this estimate is based on the quantile function of the correlated bivariate normal distribution re-written in polar coordinates (radius and angle). See \code{\link[shotGroups]{Hoyt}}. For \code{accuracy=TRUE} or three-dimensional data, it is based on the (offset) circle/sphere probabilities for the correlated multivariate normal distribution (DiDonato & Jarnagin, 1961; DiDonato, 1981). See \code{\link[shotGroups]{pmvnEll}}.
\item \code{GrubbsPearson}: The Grubbs-Pearson estimate (Grubbs, 1964) is based on the Pearson three-moment central chi^2-approximation of the true cumulative distribution function of radial error. The eigenvalues of the covariance matrix of shot-coordinates are used as variance estimates since they are the variances of the principal components (the PCA-rotated = decorrelated data).
\item \code{GrubbsPatnaik}: The Grubbs-Patnaik estimate (Grubbs, 1964) differs from the Grubbs-Pearson estimate insofar as it is based on the Patnaik two-moment central chi^2-approximation of the true cumulative distribution function of radial error.
\item \code{GrubbsLiu}: The Grubbs-Liu estimate was not proposed by Grubbs but follows the same principle as his original estimates. It differs from them insofar as it is based on the Liu-Tang-Zhang four-moment non-central chi^2-approximation of the true cumulative distribution function of radial error. For \code{accuracy=FALSE}, it is identical to \code{GrubbsPearson}.
\item \code{Rayleigh}: This estimate uses the Rayleigh quantile function with a bias-corrected estimate of its scale parameter (Singh, 1992). The Rayleigh distribution is valid for uncorrelated multivariate normal coordinates with equal variances.
}
}
\value{
A vector with the (chosen) hit-probability estimates.
}
\references{
DiDonato, A. R. (1988). Integration of the trivariate normal distribution over an offset spehere and an inverse problem (NSWC TR 87-27). Dahlgren, VA: U.S. Naval Surface Weapons Center Dahlgren Division.
  
DiDonato, A. R., & Jarnagin, M. P. (1961). Integration of the general bivariate Gaussian distribution over an offset circle. Mathematics of Computation, 15 (76), 375-382.

Grubbs, F. E. (1964). Approximate circular and noncircular offset probabilities of hitting. Operations Research, 12(1), 51-62.

Liu, H., Tang, Y., & Zhang, H. H. (2009). A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables. Computational Statistics & Data Analysis, 53(4), 853-856.

Singh, H. P. 1992. Estimation of Circular Probable Error. The Indian Journal of Statistics, Series B 5(3), 289-305.
}
\seealso{
\code{\link[shotGroups]{Hoyt}},
\code{\link[shotGroups]{mvnEll}},
\code{\link[shotGroups]{getHoytParam}},
\code{\link[shotGroups]{getRayParam}},
\code{\link[shotGroups]{getCEP}},
\code{\link[robustbase]{covMcd}},
\code{\link[shotGroups]{getConfEll}}
}
\examples{
# coordinates given by a suitable data frame
# estimated fraction of shots within a circle with radius 1 MOA.
getHitProb(DFscar17, r=1, unit='MOA', accuracy=FALSE,
           dstTarget=100, conversion='yd2in',
           type=c('CorrNormal', 'GrubbsPatnaik'))

# coordinates given by a matrix
\dontrun{
xy <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
getHitProb(xy, r=1, unit='MOA', accuracy=FALSE,
           dstTarget=100, conversion='yd2in',
           type=c('CorrNormal', 'GrubbsPatnaik'))
}
}
