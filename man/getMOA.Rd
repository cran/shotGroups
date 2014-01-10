\name{getMOA}
\alias{getMOA}
\title{Conversion to minutes of angle}
\description{Converts object size to MOA (minutes of angle).}
\usage{getMOA(x, dst, conversion = 'm2cm')}
\arguments{
  \item{x}{a numerical vector of object sizes.}
  \item{dst}{a numerical vector of viewing distances.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates. Either a character vector indicating the conversion like \code{'m2cm'} for meter to centimeter, \code{'yd2in'} for yards to inches, or \code{'ft2cm'} for feet to cm. Imperial (yd, ft, in) and metric units (m, cm, mm) can be freely mixed. Alternatively, a numerical vector giving the multiplication factor for conversion: 100 for m to cm, 36 for yd to in, and 12 for ft to in.}
}
\value{A numerical vector with the MOA values.}
\seealso{\code{\link[shotGroups]{fromMOA}}}
\examples{
objSize <- seq(10, 30, by=5)  # cm
objDst  <- 25                 # m
getMOA(objSize, dst=objDst, conversion='m2cm')

# this should return objSize
fromMOA(getMOA(objSize, objDst, conversion='m2cm'),
        objDst, conversion='m2cm')
}
