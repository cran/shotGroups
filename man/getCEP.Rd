\name{getCEP}
\alias{getCEP}
\alias{getCEP.data.frame}
\alias{getCEP.default}
\title{Circular error probable (CEP)}
\description{Calculates the Rayleigh, Grubbs-Patnaik, and modified RAND R-234 estimate for the circular error probable. CEP is defined as the radius of the circle to contain a certain fraction of points of impact, typically 50\%.}
\usage{getCEP(xy, dstTarget = 100, conversion = 'm2cm', accuracy = FALSE,
       type = c('Rayleigh', 'Grubbs', 'RAND'))

\method{getCEP}{data.frame}(xy, dstTarget = 100, conversion = 'm2cm', accuracy = FALSE,
       type = c('Rayleigh', 'Grubbs', 'RAND'))

\method{getCEP}{default}(xy, dstTarget = 100, conversion = 'm2cm', accuracy = FALSE,
       type = c('Rayleigh', 'Grubbs', 'RAND'))}
\arguments{
  \item{xy}{either a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point), or a data frame with either the variables \code{X}, \code{Y} or \code{Point.X}, \code{Point.Y}.}
  \item{dstTarget}{a numerical value with the distance to the target - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates - used in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{accuracy}{a logical value indicating whether systematic location bias should be taken into account (see details).}
  \item{type}{string vector indicating which CEP estimate to report.}
}
\details{No coordinate transforms are done (unlike in \code{\link[shotGroups]{groupLocation}}, \code{\link[shotGroups]{groupShape}}, \code{\link[shotGroups]{groupSpread}}), i.e., origin is not assumed to be in top-left corner, and points are not taken relative to point of aim. \cr \cr
By default, the reported CEP estimates do not take into account accuracy, i.e., any systematic location bias. They are then based on the centered data, and are thus a pure precision (spread) measure. Set \code{accuracy=TRUE} to incorporate systematic deviation from point of aim.
\itemize{
\item The Rayleigh 50\%-CEP estimate uses the quantile function of the Rayleigh distribution with a bias-corrected estimate of its scale parameter (Singh, 1992).
\item The Grubbs-Patnaik 50\%-CEP estimate is based on the quantile function of a chi-square distribution (Williams, 1997). The eigenvalues of the covariance matrix of (x,y)-coordinates are used here as they are the variances of the principal components (the PCA-rotated = decorrelated data).
\item The modified RAND R-234 estimate for the 50\%-CEP is a weighted sum of the square root of the aforementioned eigenvalues, that is of the standard deviations of the two principal components. The bias correction is from Williams (1997). For zero correlation and equal variances of x- and y-coordinates, the modified RAND estimate is very similar to the radius of the parametric 50\%-confidence ellipse under the assumption of bivariate normality.
}
The estimates for 90\% and 95\% are based on linear interpolation of 50\% values according to table 3 in McMillan and McMillan (2008) to correct for the aspect ratio of the error ellipse. Warnings are given if the data shows ellipticity or bias greater than what CEP is deemed useful for.
}
\value{
A list with the (chosen) CEP estimates.
\item{Rayleigh}{a numerical matrix with the Rayleigh estimate for the CEP at levels 50\%, 90\%, and 95\% (in original measurement units and MOA).}
\item{Grubbs}{a numerical matrix with the Grubbs-Patnaik estimate for the CEP at levels 50\%, 90\%, and 95\% (in original measurement units and MOA).}
\item{RAND}{a numerical matrix with the modified RAND R-234 estimate for the CEP at levels 50\%, 90\%, and 95\% (in original measurement units and MOA).}
\item{ellShape}{a numerical vector with the aspect ratio and the flattening of the error ellipse.}
\item{ctr}{group center}
}
\references{
McMillan, C. and McMillan, P. 2008. Characterizing rifle performance using circular error probable measured via a flatbed scanner. \url{http://statshooting.com/}

Singh, H. P. 1992. Estimation of Circular Probable Error. The Indian Journal of Statistics, Series B 54(3), 289-305.

Williams, C. E. 1997. A comparison of circular error probable estimators for small samples. M.Sc. thesis. AFIT/GOA/ENS/97M-14.
}
\seealso{
\code{\link[shotGroups]{getRayParam}},
\code{\link[shotGroups]{getConfEll}}
}
\examples{
# coordinates given by a suitable data frame
(cep <- getCEP(DFtalon, accuracy=FALSE, dstTarget=10, conversion='m2mm',
               type=c('Rayleigh', 'Grubbs', 'RAND')))

# plot points, centers, and circles indicating 50\%-CEP estimates
plot(Point.Y ~ Point.X, data=DFtalon, asp=1, pch=16)
drawCircle(cep$ctr, cep$Rayleigh['unit', '50\%'], fg='red')
drawCircle(cep$ctr, cep$Grubbs['unit',   '50\%'], fg='green3')
drawCircle(cep$ctr, cep$RAND['unit',     '50\%'], fg='blue')
points(cep$ctr[1], cep$ctr[2], pch=4, col='gray50', cex=2, lwd=2)
legend(x='bottomleft',
       legend=c('Rayleigh 50', 'Grubbs 50', 'RAND 50', 'center'),
       col=c('red', 'green3', 'blue', 'gray50'), lwd=2,
       lty=c(1, 1, 1, NA), pch=c(NA, NA, NA, 4), bg='white')

# calculate actual coverage percentage of 50% CEP estimates
dists <- getDistToCtr(DFtalon)           # distances to center

# extract CEP radius for all estimates
CEPr  <- c(Rayleigh=cep$Rayleigh['unit', '50\%'],
             Grubbs=cep$Grubbs['unit',   '50\%'],
               RAND=cep$RAND['unit',     '50\%'])

# percentage of points in circle with radius = CEP
100 * sapply(CEPr, function(x) sum(dists <= x)) / length(dists)

# coordinates given by a matrix
\dontrun{
xy <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
getCEP(xy, accuracy=FALSE, dstTarget=25, conversion='m2cm',
       type=c('Rayleigh', 'Grubbs', 'RAND')))
}
}
