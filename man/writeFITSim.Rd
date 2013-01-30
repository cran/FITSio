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
crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA,
axDat = NA, header = '')

writeFITSim16i(X, file = "R.fits", ...)
}

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
  \item{axDat}{Data frame with axis data (optional), see details.}
  \item{header}{Header (optional) as 80-character card images.}
  \item{\dots}{Arguments as defined above, as needed, for
          \code{writeFITSim16i}.}
}
\details{
  \code{writeFITSim} and \code{writeFITSim16i} write multi-dimensional
  data arrays and
  header information to FITS-format files.  A single image is a
  two-dimensional array; data cubes contain two dimensions plus one
  additional dimension for each (often velocity) plane.

  Axis data may be given as a data frame \code{axDat}, typically from
  \code{\link{readFITS}} or, if that is missing, by individual vectors for
  \code{crpixn, crvaln}, etc.  \code{axDat} has priority over individual
  vectors.  \code{axDat} may be edited easily with e.g.
  \code{axDat <- edit(axDat)}.

  Header data may be added to that the function automatically generates.
  No or little editing is needed if the header is taken from an existing
  file by e.g. \code{\link{readFITSheader}}: \code{writeFITSim} will
  remove reserved keywords in the header that conflict with those
  generated and add the END statement.  See \code{\link{modifyHeader}}
  for functions to add to or modify the header.

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

\seealso{\code{\link{readFITS}}, \code{\link{readFITSheader}},
     \code{\link{modifyHeader}}}

\examples{

## Make data array with axis information, write to disk
X <- matrix(1:15, ncol = 3)
writeFITSim(X, file = "test.fits", c1 = "Test FITS file",
            crpixn = c(1,1), crvaln = c(10, 100), cdeltn = c(8, 2),
            ctypen = c("Distance", "Time"),
            cunitn = c("Furlongs", "Fortnights"))
## Read back in, modify data offset, header, and axis data information,
## then write modified version as new file
Z <-  readFITS("test.fits")
Z$imDat <- Z$imDat + 300
Z$header <- addKwv('SCALE', 1.03, 'test header mod', header=Z$header)
# Z$axDat <- edit(Z$axDat)  # interactive edits
Z$axDat$cdelt
Z$axDat$cdelt[2] <- 20
Z$axDat$cdelt
writeFITSim(Z$imDat, file='test_modified.fits', axDat=Z$axDat, header=Z$header)
unlink("test.fits")

### 3-dimensional array example
## Write sample file
X <- array(1:210, dim = c(10, 7, 3))
writeFITSim(X, "test.fits")
## Read back in and display plane 2, no axis scale markings
Z <-  readFITS("test.fits")
dim(Z$imDat[,,2])
image(Z$imDat[,,2], xaxt = "n", yaxt = "n")
str(Z)
print(Z$axDat)
unlink("test.fits")

### Note: either of the writeFITSim() calls here could be replaced
###       with writeFITSim16i() calls with the identical argument.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
