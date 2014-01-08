\name{readDataMisc}
\alias{readDataMisc}
\title{Read data from text files}
\description{Reads data from text files that have a similar structure to OnTarget PC/TDS output files. Several files can be read with one call.}
\usage{readDataMisc(fPath = getwd(), fNames, fPat, combine = TRUE)}
\arguments{
  \item{fPath}{a character string containing the path to the folder with the data files, e.g. \code{'c:/folder/otFiles'}.}
  \item{fNames}{a character vector containing the file names of the files that should be read in.}
  \item{fPat}{a character string containing the regular-expression that describes all names of files that should be read in. E.g., \code{'^points[[:digit:]]{2}\\\\.txt$'} for filenames \code{'points**.txt'}, where \code{**} are 2 digits. See \code{\link{regex}}, \code{\link[utils]{glob2rx}}.}
  \item{combine}{a logical value indicating whether the data should be combined into one big data frame with \code{\link[shotGroups]{combineData}}.}
}
\details{
If \code{fNames} is provided, \code{fPat} is ignored. \cr \cr
If neither \code{fNames} nor \code{fPat} is provided, and we are in interactive mode under Windows, files can be chosen interactively. \cr \cr
This function is basically a wrapper for \code{\link[utils]{read.table}}. \cr \cr
The files need to be whitespace-delimited, and have a header with the variable names. Variable names must not contain spaces. They need to contain at least the variable \code{Group} if \code{\link[shotGroups]{combineData}} should be later used to combine them into one big data frame. In order to be later used by functions \code{\link[shotGroups]{analyzeGroup}} or \code{\link[shotGroups]{compareGroups}}, coordinates for points of impact must be given either with \code{Point.X}, \code{Point.Y} or \code{X}, \code{Y}. Point of aim can be given with \code{Aim.X}, \code{Aim.Y}, otherwise (0,0) will be assumed in analysis functions. Distance to target can be given with \code{Distance}, otherwise 100m will be assumed in analysis functions. Additional variable names must not contain spaces. There can only be exactly as many variable names as there are non-empty data-columns. \cr \cr
For reading in files exported from OnTarget PC 1.1*, see \code{\link[shotGroups]{readDataOT1}}. \cr
For reading in files exported from OnTarget PC 2.1* or OnTarget TDS v3.7*, see \code{\link[shotGroups]{readDataOT2}}.
}
\value{
With \code{combine=FALSE}: a list of data frames, each from one file that was read in. This list can then be combined into one big data frame by \code{\link[shotGroups]{combineData}}.
\item{file1}{data frame containing data from the first file}
\item{file2}{data frame containing data from the second file}
\item{\ldots}{more data frames}
With \code{combine=TRUE}: a combined data frame from \code{\link[shotGroups]{combineData}}.
}
\seealso{
\code{\link[utils]{read.table}},
\code{\link{regex}},
\code{\link[utils]{glob2rx}},
\code{\link[shotGroups]{combineData}},
\code{\link[shotGroups]{readDataOT1}},
\code{\link[shotGroups]{readDataOT2}}
}
\examples{
\dontrun{
fPath   <- 'c:/folder/files'             # folder with data files
fileNam <- c('pts01.txt', 'pts02.txt')   # desired files in that folder
DFgroup <- readDataMisc(fPath, fNames=fileNam, combine=TRUE)

## alternatively, specify filename pattern for all files to be read in
fPath   <- 'c:/folder/otFiles'           # folder with data files
fPat    <- '^pts[[:digit:]]{2}\\\\.txt$' # filename pattern
DFgroup <- readDataMisc(fPath, fPat=pattern, combine=TRUE)
}

## result should look like this
data(DFcm)
head(DFcm)
}
