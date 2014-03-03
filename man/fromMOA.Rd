\name{fromMOA}
\alias{fromMOA}
\title{Conversion from MOA, SMOA or mil to absolute size}
\description{Converts angular measure (minute of angle MOA = arcminute, Shooter's MOA SMOA, milliradian milrad) to object size.}
\usage{fromMOA(x, dst, conversion = 'm2cm', type = c('MOA', 'SMOA', 'milrad'))}
\arguments{
  \item{x}{a numerical vector of MOA values.}
  \item{dst}{a numerical vector of viewing distances.}
  \item{conversion}{how to convert the measurement unit for distance to target to that of the (x,y)-coordinates. Either a character vector indicating the conversion like \code{'m2cm'} for meter to centimeter, \code{'yd2in'} for yards to inches, or \code{'ft2cm'} for feet to cm. Imperial (yd, ft, in) and metric units (m, cm, mm) can be freely mixed. Alternatively, a numerical vector giving the multiplication factor for conversion: 100 for m to cm, 36 for yd to in, and 12 for ft to in.}
  \item{type}{type of angular measure: \code{'MOA'} for minute of angle, \code{'SMOA'} for Shooter's MOA, \code{'milrad'} for milliradian. See details.}
}
\details{
1 MOA (minute of angle, arcmin) = 1/60 degree. Shooter's MOA = SMOA = Inches Per Hundred Yards IPHY. 1 inch at 100 yards = 1 SMOA. 1 milliradian = 1/1000 radian. Details are given in the vignette, see \code{vignette('shotGroups') .}
}
\value{A numerical vector with the object sizes. The measurement unit is determined by \code{conversion}.}
\seealso{\code{\link[shotGroups]{getMOA}}}
\examples{
size <- seq(1, 20, by=5)   # inch
dst  <- 100                # yard
fromMOA(size, dst=dst, conversion='yd2in', type='MOA')

# this should return objSize
MOA <- getMOA(size, dst=dst, conversion='yd2in', type='MOA')
fromMOA(MOA, dst=dst, conversion='yd2in', type='MOA')
        
# SMOA
fromMOA(c(1, 2, 5), dst=100, conversion='yd2in', type='SMOA')

# milliradian
fromMOA(c(1, 2, 5), dst=100, conversion='m2mm', type='milrad')
}
