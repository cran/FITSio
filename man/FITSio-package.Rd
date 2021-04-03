\name{FITSio-package}
\alias{FITSio-package}
\alias{FITSio}
\docType{package}
\title{FITS file input-output functions
}
\description{
FITS, the Flexible
Image Transport System, is
a standard data format for astronomy.
This package contains
functions to read files with FITS file image and binary and ASCII table formats
and to write FITS image files.  The functions comply with the standards
defined by the International Astronomical Union's FITS Working Group.
}
\details{
\tabular{ll}{
Package: \tab FITSio\cr
Type: \tab Package\cr
Version: \tab 2.1-6\cr
Date: \tab 2021-03-28\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
At the high level, this package contains functions to read
single FITS Header and Data Units (HDUs) containing image and binary
table extensions, and one to write FITS image
files.  \code{readFITS} automatically
recognizes
image (multi-dimensional arrays) and binary table extensions, returning a
list with data, header, and scaling information.
\code{readFrameFromFITS} returns an \R data frame from a single binary
table HDU.  Both functions accept an argument to pick out the
\emph{n}th HDU in a larger file.
Binary table complex, and array
descriptor data types are not implemented in this release due to a lack
of examples for testing.  8, 16, 24, and 32-bit bit arrays return as integers.  Other lengths
  are not supported.

For writing, the package contains the high-level function
\code{writeFITSim} to write FITS
images, along with a minor extension to efficiently write
data in 16-bit integers, \code{writeFITSim16i}.

A set of mid-level functions enables reading files with
combinations of HDUs: \code{readFITSheader} reads headers,
\code{readFITSarray} reads image extensions, and
\code{readFITSbintable} reads binary table extensions.
\code{readFITStable} reads ASCII table extensions.
\code{readFITS} and \code{readFrameFromFITS} invoke these functions to
do the actual reading.

Function \code{axVec} generates a vector for image axis labeling from
data in an image FITS file.
}
\author{Andrew Harris <harris@astro.umd.edu>, with contributions from
Eric H. Neilsen, Jr. and Bi Cheng Wu
}
\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{https://fits.gsfc.nasa.gov/}
}
\examples{
require(FITSio)
## Make test image with axis information, write to disk
Z <- matrix(1:15, ncol = 3)
filename <- paste(tempdir(), "test.fits", sep="")
writeFITSim(Z, file = filename, c1 = "Test FITS file",
            crpix = c(1,1), crvaln = c(10, 100), cdeltn = c(8, 2),
            ctypen = c("Distance", "Time"),
            cunitn = c("Furlongs", "Fortnights"))

## Read image back and display
X <-  readFITS(filename)
ax1 <- axVec(1, X$axDat)          # Make axis vector for image
ax2 <- axVec(2, X$axDat)
xlab <- X$axDat$ctype[1]
ylab <- paste(X$axDat$ctype[2], " [", X$axDat$cunit[2], "]", sep = "")
image(ax1, ax2, X$imDat, xlab = xlab, ylab = ylab)
str(X)
X$axDat                           # Display data frame with axis data
X$hdr[1:10]                       # Header sample
X$hdr[which(X$hdr=="BITPIX")+1]   # BITPIX value from header

### Read back in, modify data, header, and axis data information,
## then write modified version as new file
Z <-  readFITS(filename)
  Z$imDat <- Z$imDat + 300
  Z$header <- addKwv('SCALE', 1.03, 'test header mod', header=Z$header)
  # Z$axDat <- edit(Z$axDat)  # interactive edit
  Z$axDat$cdelt[2] <- 20
filename2 <- paste(tempdir(), "test.fits", sep="")
writeFITSim(Z$imDat, file=filename2, axDat=Z$axDat, header=Z$header)

## Clean up files to avoid clutter
unlink(filename)
unlink(filename2)
}

\keyword{file}
