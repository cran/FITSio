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
  \item{zz}{File handle; see Example.}
  \item{hdr}{Header card images, raw or parsed.}
}
\details{
  \code{readFITSarray} reads the data from the image part of a FITS
  Header and Data Unit (hdu) containing image data.  The header must be
  read first by \code{\link{readFITS}} or \code{\link{readFITSheader}};
  either this header or the
  parsed version from \code{\link{parseHdr}} are valid for the \var{hdr}
  variable.
}
\value{
  A list containing
  \item{imDat}{Data array.}
  \item{axDat}{Data frame with axis scaling and labels.}
  \item{hdr}{Vector with the parsed header.}
}

\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{https://fits.gsfc.nasa.gov/}
}

\author{Andrew Harris}

\note{
Function assigns values of 1 to CRPIX, CRVAL, and CDELT if they are
unspecified in the header.

Graphical FITS viewers such as \emph{fv}
  (\url{https://heasarc.gsfc.nasa.gov/ftools/fv/}) and \emph{SAOImage DS9}
(\url{http://ds9.si.edu/})
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
    require(FITSio)
    ## Make a test file.
    Z <- matrix(1:15, ncol = 3)
    filename <- paste(tempdir(), "test.fits", sep="")
    writeFITSim(Z, filename)
    ## Open file, read header and array, close file and delete.
    zz <- file(description = filename, open = "rb")
    header <- readFITSheader(zz)  # image data off primary header
    D <- readFITSarray(zz, header)
    close(zz)
    ## Look at data list, header file, and parsed header
    str(D)
    image(D$imDat)
    str(header)
    str(parseHdr(header))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
