\name{getCEP}
\alias{getCEP}
\alias{getCEP.data.frame}
\alias{getCEP.default}
\title{Circular Error Probable (CEP) and Spherical Error Probable (SEP)}
\description{Estimates the Circular Error Probable (CEP) or the Spherical Error Probable (SEP). CEP/SEP is defined as the radius of the circle/sphere around the point of aim such that it contains a certain fraction of points of impact, e.g., 50\% or 90\%.}
\usage{getCEP(xy, CEPlevel=0.5, dstTarget = 100, conversion = 'm2cm',
       accuracy = FALSE, type = 'CorrNormal', doRob = FALSE)

\method{getCEP}{data.frame}(xy, CEPlevel=0.5, dstTarget = 100, conversion = 'm2cm',
       accuracy = FALSE, type = 'CorrNormal', doRob = FALSE)

\method{getCEP}{default}(xy, CEPlevel=0.5, dstTarget = 100, conversion = 'm2cm',
       accuracy = FALSE, type = 'CorrNormal', doRob = FALSE)}
\arguments{
  \item{xy}{either a numerical (n x p)-matrix with the coordinates of n points (1 row of p coordinates per point), or a data frame with either the variables \code{X}, \code{Y} or \code{Point.X}, \code{Point.Y} (optionally \code{Z} or \code{Point.Z}).}
  \item{CEPlevel}{a numerical value with the coverage for the CEP/SEP.}
  \item{dstTarget}{a numerical value with the distance to the target - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the shot-coordinates - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{accuracy}{logical: take systematic location bias into account? See details.}
  \item{type}{string vector indicating which CEP/SEP estimate to report. Possible values are \code{'CorrNormal', 'GrubbsPearson', 'GrubbsPatnaik', 'GrubbsLiu', 'Rayleigh', 'Ethridge', 'RAND'}. See details.}
  \item{doRob}{logical: use robust estimation of center and covariance matrix as basis for estimators?}
}
\details{
For \code{accuracy=FALSE} (default), the reported CEP/SEP estimates do not take into account accuracy, i.e., any systematic location bias. The data is then first centered on the empirical group mean, assumed to coincide with the point of aim. The resulting CEP/SEP-MPI (around Mean Point of Impact) is a pure precision (spread) measure. Set \code{accuracy=TRUE} to incorporate systematic accuracy bias such that the point of aim is in the origin 0, possibly offset from the true group center.
\itemize{
\item \code{CorrNormal}: For \code{accuracy=FALSE} and two-dimensional data, this estimate is based on the correlated bivariate normal distribution re-written in polar coordinates (radius and angle) (see \code{\link[shotGroups]{Hoyt}}). For \code{accuracy=TRUE} or three-dimensional data, it is based on the (offset) circle/sphere probabilities for the correlated multivariate normal distribution (DiDonato & Jarnagin, 1961; DiDonato, 1981, see \code{\link[shotGroups]{qmvnEll}}). This estimate is available for all probability levels.
\item \code{GrubbsPearson}: The Grubbs-Pearson estimate (Grubbs, 1964) is based on the Pearson three-moment central chi^2-approximation of the true cumulative distribution function of radial error. The eigenvalues of the covariance matrix of shot-coordinates are used as variance estimates since they are the variances of the principal components (the PCA-rotated = decorrelated data). This estimate is available for all probability levels, and generalizes to three dimensions.
\item \code{GrubbsPatnaik}: The Grubbs-Patnaik estimate (Grubbs, 1964) differs from the Grubbs-Pearson estimate insofar as it is based on the Patnaik two-moment central chi^2-approximation of the true cumulative distribution function of radial error.
\item \code{GrubbsLiu}: The Grubbs-Liu estimate was not proposed by Grubbs but follows the same principle as his original estimates. It differs from them insofar as it is based on the Liu-Tang-Zhang four-moment non-central chi^2-approximation of the true cumulative distribution function of radial error. For \code{accuracy=FALSE}, it is identical to \code{GrubbsPearson}.
\item \code{Rayleigh}: For \code{accuracy=FALSE} and two-dimensional data, this estimate uses the Rayleigh distribution (see \code{\link[shotGroups]{getRayParam}}). It is valid for uncorrelated bivariate normal coordinates with equal variances. This estimate is available for all probability levels. For \code{accuracy=FALSE} and three-dimensional data, the Maxwell-Boltzmann distribution is used (see \code{\link[shotGroups]{getMaxParam}}). For \code{accuracy=TRUE} and two-dimensional data, the estimate uses the Rice distribution (see \code{\link[shotGroups]{getRiceParam}}). For \code{accuracy=TRUE} and three-dimensional data, it is based on the offset sphere probabilities for the multivariate normal distribution set to have equal variances (see \code{\link[shotGroups]{qmvnEll}}).
\item \code{Ethridge}: The Ethridge estimate (Ethridge, 1983; Puhek, 1992) is not based on the assumption of multivariate normality of coordinates but uses a robust unbiased estimator for the median radius (Hogg, 1967). It can only be reported for probability 0.5 but generalizes to three dimensions.
\item \code{RAND}: The modified RAND R-234 estimate (RAND, 1952; Pesapane & Irvine, 1977; Puhek 1992) is a weighted sum of the square root of the eigenvalues of the covariance matrix, that is of the standard deviations of the two principal components. It can only be reported for probability 0.5 and does not generalize to three dimensions.
}
}
\value{
A list with the (chosen) CEP/SEP estimates and supplementary information on the distribution shape.
\item{CEP}{a numerical matrix with the chosen CEP/SEP estimates at the indicated \code{CEPlevel} (in original measurement units, MOA, SMOA, milliradian).}
\item{ellShape}{a numerical vector with the aspect ratio of the ellipse (the square root of condition index \code{\link{kappa}}) and its flattening (1 - inverse aspect ratio).}
\item{ctr}{group center}
}
\references{
DiDonato, A. R. (1988). Integration of the trivariate normal distribution over an offset spehere and an inverse problem (NSWC TR 87-27). Dahlgren, VA: U.S. Naval Surface Weapons Center Dahlgren Division.

DiDonato, A. R., & Jarnagin, M. P. (1961). Integration of the general bivariate Gaussian distribution over an offset circle. Mathematics of Computation, 15 (76), 375-382.

Grubbs, F. E. (1964). Approximate circular and noncircular offset probabilities of hitting. Operations Research, 12(1), 51-62.

Hogg, R. V. (1967). Some observations on robust estimation. Journal of the American Statistical Association, 62 (320), 1179-1186.

Liu, H., Tang, Y., & Zhang, H. H. (2009). A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables. Computational Statistics & Data Analysis, 53(4), 853-856.

Pesapane, J., & Irvine, R. B. (1977). Derivation of CEP formula to approximate RAND-234 tables. Offut AFB, NE: Ballistic Missile Evaluation, HQ SAC.

Puhek, P. (1992). Sensitivity analysis of circular error probable approximation techniques (AFIT/GOR/ENS/92M-23). Wright-Patterson AFB, OH: Air Force Institute of Technology.

RAND Corporation. (1952). Offset circle probabilities (RAND-234). Santa Monica, CA: RAND Corporation.

Singh, H. P. 1992. Estimation of Circular Probable Error. The Indian Journal of Statistics, Series B 5(3), 289-305.
}
\seealso{
\code{\link[shotGroups]{Rayleigh}},
\code{\link[shotGroups]{Maxwell}},
\code{\link[shotGroups]{Hoyt}},
\code{\link[shotGroups]{Rice}},
\code{\link[shotGroups]{mvnEll}},
\code{\link[shotGroups]{getHoytParam}},
\code{\link[shotGroups]{getRayParam}},
\code{\link[shotGroups]{getRiceParam}},
\code{\link[shotGroups]{getConfEll}},
\code{\link[shotGroups]{getHitProb}},
\code{\link[robustbase]{covMcd}}
}
\examples{
# coordinates given by a suitable data frame
(cep <- getCEP(DFtalon, CEPlevel=0.5, accuracy=FALSE,
               dstTarget=10, conversion='m2mm',
               type=c('CorrNormal', 'GrubbsPatnaik', 'Rayleigh')))

# plot points, centers, and circles indicating 50\%-CEP estimates
plot(point.y ~ point.x, data=DFtalon, asp=1, pch=16)
drawCircle(cep$ctr, cep$CEP['unit',    'CorrNormal'], fg='red')
drawCircle(cep$ctr, cep$CEP['unit', 'GrubbsPatnaik'], fg='green3')
drawCircle(cep$ctr, cep$CEP['unit',      'Rayleigh'], fg='blue')
points(cep$ctr[1], cep$ctr[2], pch=4, col='gray50', cex=2, lwd=2)
legend(x='bottomleft',
       legend=c('Grubbs 50', 'Corr Normal 50', 'Rayleigh 50', 'center'),
       col=c('red', 'green3', 'blue', 'gray50'), lwd=2,
       lty=c(1, 1, 1, NA), pch=c(NA, NA, NA, 4), bg='white')

# calculate actual coverage percentage of 50% CEP estimates
dists <- getDistToCtr(DFtalon)           # distances to center

# extract CEP radius for all estimates
CEPr <- c(Grubbs=cep$CEP['unit', 'GrubbsPatnaik'],
      CorrNormal=cep$CEP['unit',    'CorrNormal'],
        Rayleigh=cep$CEP['unit',      'Rayleigh'])

# percentage of points in circle with radius = CEP
100 * sapply(CEPr, function(x) sum(dists <= x)) / length(dists)

# coordinates given by a matrix
\dontrun{
xy <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
getCEP(xy, accuracy=FALSE, dstTarget=25, conversion='m2cm',
       type=c('Grubbs', 'CorrNormal', 'Rayleigh')))
}
}
