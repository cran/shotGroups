\name{getCEP}
\alias{getCEP}
\title{Circular error probable (CEP)}
\description{Gives the Rayleigh and the modified RAND R-234 estimate for the circular error probable. CEP is defined as the radius of the circle to contain a certain fraction of points of impact.}
\usage{getCEP(xy, dstTarget = 25, conversion = 'm2cm')}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point).}
  \item{dstTarget }{a numerical value with the distance to the target in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates in MOA calculation. See \code{\link[shotGroups]{getMOA}}.}
}
\details{The 50\%-CEP estimate based on the Rayleigh distribution is a constant fraction of the mean distance to center. The modified RAND R-234 estimate for the 50\%-CEP is a weighted sum of the square root of the eigenvalues of the covariance matrix of (x,y)-coordinates. These are the standard deviations of the two principal components of the data. For zero correlation and equal variances of x- and y-coordinates, the RAND estimate is very similar to the radius of the parametric 50\%-confidence ellipse under the assumption of bivariate normality. The estimates for 90\% and 95\% are based on table 3 in the reference given below.}
\value{
A list with the CEP estimates.
\item{RAND }{a numerical matrix with the modified RAND R-234 estimate for the CEP at levels 50\%, 90\%, and 95\% (in original measurement units and MOA).}
\item{Rayleigh }{a numerical matrix with the Rayleigh estimate for the CEP at levels 50\%, 90\%, and 95\% (in original measurement units and MOA).}
\item{ellShape }{a numerical vector with the aspect ratio and the flattening of the error ellipse.}
}
\references{
McMillan, C. and McMillan, P. 2008. Characterizing rifle performance using circular error probable measured via a flatbed scanner. \url{http://statshooting.com/}
}
\note{The Rayleigh approximation emerged under historical circumstances for ease of computation. The regular confidence ellipse as returned by \code{\link[shotGroups]{getConfEll}} should be considered instead.}
\seealso{\code{\link[shotGroups]{getConfEll}}}
\examples{
xy <- matrix(round(rnorm(100, 0, 5), 2), ncol=2)
(cep <- getCEP(xy, dstTarget=25, conversion='m2cm'))

# plot points and circle indicating 90\%-CEP
ctr <- colMeans(xy)             # center point
plot(xy, pch=16, asp=1)
drawCircle(ctr[1], ctr[2], cep$RAND["unit", "90\%"], fg='blue')
drawCircle(ctr[1], ctr[2], cep$Rayleigh["unit", "90\%"], fg='green3')
points(ctr[1], ctr[2], pch=4, col="red", cex=2)
legend(x="bottomleft", legend=c("RAND 90", "Rayleigh 90", "center"),
       col=c("blue", "green3", "red"), lwd=c(2, 2, NA),
       lty=c(1, 1, NA), pch=c(NA, NA, 4))
}
