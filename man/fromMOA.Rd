\name{fromMOA}
\alias{fromMOA}
\title{Conversion from minutes of angle}
\description{Converts MOA (minutes of angle) to object size.}
\usage{fromMOA(x, dst, conversion = 'm2cm')}
\arguments{
  \item{x }{a numerical vector of MOA values.}
  \item{dst }{a numerical vector of viewing distances.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates. Either a character vector indicating the conversion: \code{'m2cm'} for meter to centimeter, \code{'m2mm'} for meter to millimeter, \code{'yd2in'} for yards to inches, or \code{'ft2in'} for feet to inches. Alternatively, a numerical vector giving the multiplication factor for conversion: 100 for m to cm, 36 for yd to in, and 12 for ft to in.}
}
\value{A numerical vector with the object sizes. The measurement unit is determined by \code{conversion}.}
\section{Warning }{The calculation leads to rounding errors after a few significant digits.}
\seealso{\code{\link[shotGroups]{getMOA}}}
\examples{
objMOA <- seq(10, 30, by=5)
objDst <- 25                 # m
fromMOA(objMOA, dst=objDst, conversion='m2cm')

## check rounding problems
objSize <- seq(10, 30, by=5)  # cm

## this should return objSize
fromMOA(getMOA(objSize, objDst, conversion='m2cm'),
        objDst, conversion='m2cm')
}
