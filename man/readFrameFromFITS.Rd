\name{readFrameFromFITS}
\alias{readFrameFromFITS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a single data set from a FITS file into a data frame}
\description{
  Read a binary table
  from a FITS file directly into an \R data frame.  
}
\usage{
readFrameFromFITS(file, hdu = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{existing FITS file name.}
  \item{hdu}{position of Header and Data Unit (HDU) in FITS file: 1 for
    the first HDU, 2 for the second HDU, etc.}
}
\details{
  \code{readFrameFromFITS} reads a single binary table Header and Data Unit
  (HDU) pair from a file and returns the values as a data table.

  Binary table bit, complex, and array
  descriptor data types are not implemented in this release due to a lack
  of examples for testing.
}
\value{
  An \R data frame with the contents of the
  requested binary table.
}

\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{http://fits.gsfc.nasa.gov/}
}

\author{Eric H.\ Neilsen, Jr.}


\note{
  Graphical FITS viewers such as \emph{fv}
  (\url{http://heasarc.gsfc.nasa.gov/ftools/fv/}) and \emph{SAOImage DS9}
(\url{http://hea-www.harvard.edu/RD/ds9/}) 
have excellent facilities for displaying FITS data, headers, and file
structure.  Having one or more graphical viewers available will
prove extremely
useful
for working with FITS
files, even when the data are read into \R for further processing.
\emph{fv} and \emph{SAOImage DS9} are in active devlopement with support for
unix, Windows, and Mac OS-X operating systems, and are available at no cost.
}

\seealso{\code{\link{readFITS}}}
\examples{
require("FITSio")

## Either download example file from
## <http://fits.gsfc.nasa.gov/fits_samples.html>
## and use
\dontrun{filename <- "IUElwp25637mxlo.fits"}
## or, for local example use
filename <- system.file("fitsExamples", "IUElwp25637mxlo.fits",
                        package = "FITSio")

## Get data and display
F <- readFrameFromFITS(filename) 
names(F)
plot(F$NET, ylab = "Value", main = names(F)[5], type = "l")

### Simple flat file example
filename <- system.file("fitsExamples", "2008_03_11_203150.fits",
                        package = "FITSio")  
F <- readFrameFromFITS(filename) 
names(F)
attach(F)
plot(DMJD, TiltX, xlab = "Time [DMJD]", ylab = "X Tilt [degr]", type = "l")
detach(F)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
