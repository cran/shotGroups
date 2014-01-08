\name{fromMOA}
\alias{fromMOA}
\title{Conversion from minutes of angle}
\description{Converts MOA (minutes of angle) to object size.}
\usage{fromMOA(x, dst, conversion = 'm2cm')}
\arguments{
  \item{x}{a numerical vector of MOA values.}
  \item{dst}{a numerical vector of viewing distances.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates. Either a character vector indicating the conversion like \code{'m2cm'} for meter to centimeter, \code{'yd2in'} for yards to inches, or \code{'ft2cm'} for feet to cm. Imperial (yd, ft, in) and metric units (m, cm, mm) can be freely mixed. Alternatively, a numerical vector giving the multiplication factor for conversion: 100 for m to cm, 36 for yd to in, and 12 for ft to in.}
}
\value{A numerical vector with the object sizes. The measurement unit is determined by \code{conversion}.}
\seealso{\code{\link[shotGroups]{getMOA}}}
\examples{
objMOA <- seq(10, 30, by=5)
objDst <- 25                  # m
fromMOA(objMOA, dst=objDst, conversion='m2cm')

## this should return objMOA
getMOA(fromMOA(objMOA, objDst, conversion='m2cm'),
       objDst, conversion='m2cm')
}
