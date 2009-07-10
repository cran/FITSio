\name{readFITSarray}
\alias{readFITSarray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read an image (multi-dimensional array) from a FITS file}
\description{Read an image (multi-dimensional array) from an open
  connection to a FITS file.
}
\usage{
readFITSarray(zz, hdr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zz}{file handle; see Example.}
  \item{hdr}{vector of header information from \code{readFITSheader}.}
}
\details{
  \code{readFITSarray} reads the data from the image part of a FITS
  Header and Data Unit containing image data.  The header must be
  read first by \code{readFITSheader} to provide the variable \var{hdr}.
}
\value{
  A list containing
  \item{imDat}{Data array.}
  \item{axDat}{Data frame with axis scaling and labels.}
  \item{hdr}{Vector with parsed header.}
}

\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{http://fits.gsfc.nasa.gov/}
}

\author{Andrew Harris}

\note{
Function assigns values of 1 to CRPIX, CRVAL, and CDELT if they are
unspecified in the header.

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

\seealso{
  \code{\link{readFITS}}, \code{\link{readFITSheader}},
  \code{\link{readFITSbintable}}, \code{\link{file}}
}

\examples{
    require("FITSio")
    ## Make a test file.
    Z <- matrix(1:15, ncol = 3)
    writeFITSim(Z, "test.fits")
    ## Open file, read header and array, close file and delete.
    zz <- file(description = "test.fits", open = "rb")
    hdr <- readFITSheader(zz)  # image data off primary header
    D <- readFITSarray(zz, hdr)
    close(zz)
    ## Look at data
    summary(D)
    image(D$imDat)
    ## Delete test file
    unlink("test.fits")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
