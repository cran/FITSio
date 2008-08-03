\name{writeFITSim}
\alias{writeFITSim}
\alias{writeFITSim16i}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write a FITS image (multi-dimensional numeric array) to disk}
\description{
  Write a FITS image (multi-dimensional numeric array) to disk.
}
\usage{
writeFITSim(X, file = "R.fits", type = "double", 
bscale = 1, bzero = 0, c1 = NA, c2 = NA, 
crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA)

writeFITSim16i(X, file = "R.fits", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Multi-dimensional numeric data array; see Details.}
  \item{file}{Output filename.}
  \item{type}{Type to write: single or double precision.}
  \item{bscale}{Global scaling factor, FITS standard meaning.}
  \item{bzero}{Global shift, FITS standard meaning.}
  \item{c1}{Character string comment line for header.}
  \item{c2}{Character string comment line for header.}
  \item{crpixn}{Vector of reference pixel numbers for axes, FITS
    standard meaning.}
  \item{crvaln}{Vector of values at reference pixels, FITS standard meaning.}
  \item{cdeltn}{Vector of axis increments per pixel, FITS standard meaning.}
  \item{ctypen}{String vector of descriptive labels for axis, FITS
    standard meaning.}
  \item{cunitn}{String vector of physical units for axis, FITS standard
    meaning.}
  \item{\dots}{Arguments as defined above, as needed, for 
          \code{writeFITSim16i}.}
}
\details{
  \code{writeFITSim} and \code{writeFITSim16i} write multi-dimensional 
  data arrays and
  header information to FITS-format files.  A single image is a
  two-dimensional array; data cubes contain two dimensions plus one
  additional dimension for each (often velocity) plane.

  \code{writeFITSim} writes integer or float data matching the data type in
  the input array.  \code{writeFITSim16i} scales and shifts 
  the input to write
  16-bit integer data with the maximum precision allowed by word
  length.  FITS variables BSCALE and BZERO are automatically updated to
  allow reconstruction of the original data values.  In cases where full
  precision is not needed, this can reduce file sizes by a factor of
  about four compared with a double-precision float.
}

\value{
  FITS file written to disk.
}

\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{http://fits.gsfc.nasa.gov/}
}

\author{Andrew Harris}

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

## Make data array with axis information, write to disk
X <- matrix(1:15, ncol = 3)
writeFITSim(X, file = "test.fits", c1 = "Test FITS file",
            crpix = c(1,1), crvaln = c(10, 100), cdeltn = c(8, 2),
            ctypen = c("Distance", "Time"),
            cunitn = c("Furlongs", "Fortnights"))
## Read back in and display
Z <-  readFITS("test.fits")
ax1 <- axVec(1, Z$axDat)          # Make axis vector for image
ax2 <- axVec(2, Z$axDat)
xlab <- Z$axDat$ctype[1]          
ylab <- paste(Z$axDat$ctype[2], " [", Z$axDat$cunit[2], "]", sep = "")
image(ax1, ax2, Z$imDat, xlab = xlab, ylab = ylab)
summary(Z)
Z$axDat                           # Display data frame with axis data
Z$hdr[1:10]                       # Header sample
Z$hdr[which(Z$hdr=="BITPIX")+1]   # BITPIX value from header
unlink("test.fits")

### 3-dimensional array example
## Write sample file
X <- array(1:210, dim = c(10, 7, 3))
writeFITSim(X, "test.fits")
## Read back in and display plane 2, no axis scale markings
Z <-  readFITS("test.fits")
dim(Z$imDat[,,2])
image(Z$imDat[,,2], xaxt = "n", yaxt = "n")
summary(Z)
print(Z$axDat)
unlink("test.fits")

### Note: either of the writeFITSim() calls here could be replaced 
###       with writeFITSim16i() calls with the identical argument.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
