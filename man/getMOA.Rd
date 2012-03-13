\name{getMOA}
\alias{getMOA}
\title{Conversion to minutes of angle}
\description{Converts object size to MOA (minutes of angle).}
\usage{getMOA(x, dst, conversion = 'm2cm')}
\arguments{
  \item{x }{a numerical vector of object sizes.}
  \item{dst }{a numerical vector of viewing distances.}
  \item{conversion }{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates. Either a character vector indicating the conversion: \code{'m2cm'} for meter to centimeter, \code{'m2mm'} for meter to millimeter, \code{'yd2in'} for yards to inches, or \code{'ft2in'} for feet to inches. Alternatively, a numerical vector giving the multiplication factor for conversion: 100 for m to cm, 36 for yd to in, and 12 for ft to in.}
}
\value{A numerical vector with the MOA values.}
\section{Warning }{The calculation leads to rounding errors after a few significant digits.}
\seealso{\code{\link[shotGroups]{fromMOA}}}
\examples{
objSize <- seq(10, 30, by=5)  # cm
objDst  <- 25                 # m
getMOA(objSize, dst=objDst, conversion='m2cm')

## check rounding problems
## this should return objSize
fromMOA(getMOA(objSize, objDst, conversion='m2cm'),
        objDst, conversion='m2cm')
}
