\name{groupShape}
\alias{groupShape}
\title{Shape analysis for a single group of bullet holes}
\description{Assesses shape of a single group of bullet holes: Outlier analysis as well as numerical and graphical normality checks for a set of (x,y)-coordinates.}
\usage{groupShape(xy, plots = TRUE, bandW = 0.5)}
\arguments{
  \item{xy }{a numerical (n x 2)-matrix with the (x,y)-coordinates of n points (1 row of coordinates per point).}
  \item{plots }{a logical value indicating whether diagrams should be shown.}
  \item{bandW }{for argument \code{bandwith} of \code{\link[graphics]{smoothScatter}}.}
}
\details{
In addition to the numerical results listed below, this function produces the following diagrams:
\itemize{
\item a combined plot for multivariate outlier identification as produced by \code{\link[mvoutlier]{aq.plot}}
\item a chi-square Q-Q-plot for eyeballing multivariate normality as produced by \code{\link[mvoutlier]{chisq.plot}}, including a reference line with intercept 0 and slope 1
\item a heatmap of a 2D-kernel density estimate for the (x,y)-coordinates as produced by \code{\link[graphics]{smoothScatter}} together with group center and error ellipse based on a robust estimate for the covariance matrix
\item a Q-Q-plot for x-coordinates for eyeballing normality
\item a Q-Q-plot for y-coordinates for eyeballing normality
\item a histogram for x-coordinates including a fitted normal distribution as well as a nonparametric kernel density estimate
\item a histogram for y-coordinates including a fitted normal distribution as well as a nonparametric kernel density estimate
}
}
\value{
A list with the results from the numerical analyses and statistical tests.
\item{corXY }{correlation matrix of (x,y)-coordinates.}
\item{corXYrob }{robust estimate of correlation matrix of (x,y)-coordinates.}
\item{Outliers }{a vector of row indices for observations identified as outliers.}
\item{ShapiroX }{Shapiro-Wilk-Test result for normality of x-coordinates.}
\item{ShapiroY }{Shapiro-Wilk-Test result for normality of y-coordinates.}
\item{multNorm }{E-statistic-Test result for multivariate normality of (x,y)-coordinates.}
}
\seealso{
\code{\link[stats]{qqnorm}},
\code{\link[graphics]{smoothScatter}},
\code{\link[graphics]{hist}},
\code{\link[stats]{kernel}},
\code{\link[robustbase]{covMcd}},
\code{\link[stats]{shapiro.test}},
\code{\link[energy]{mvnorm.etest}},
\code{\link[mvoutlier]{chisq.plot}},
\code{\link[mvoutlier]{aq.plot}}
}
\examples{
xy  <- matrix(round(rnorm(200, 0, 5), 2), ncol=2)
res <- groupShape(xy, bandW=1.6)
names(res)
res$corXY
res$Outliers
res$multNorm
}
