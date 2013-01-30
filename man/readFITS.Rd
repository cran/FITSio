\name{readFITS}
\alias{readFITS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a single data set from a FITS file}
\description{
  Read a single image (multi-dimensional array) or binary table
  from a FITS file.  The source code \code{readFITS.r} is a model for
  creating code to read more complex FITS files.
}
\usage{
readFITS(file = "R.fits", hdu = 1, maxLines = 5000,
fixHdr = c('none', 'remove', 'substitute'), phdu = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{Existing FITS file name.}
  \item{hdu}{Position of Header and Data Unit (HDU) in FITS file: 1 for
    the first HDU, 2 for the second HDU, etc.}
  \item{maxLines}{Integer: maximum number of header lines to read.}
  \item{fixHdr}{Deal with (illegal) nonprinting characters in the
    header.}
  \item{phdu}{Rarely needed; see Details.}

}
\details{
  \code{readFITS} is a simple but complete FITS file reader,
  automatically detecting image and binary table data extensions
  (ASCII tables and random groups are
  not supported in this
  release).  It reads a single Header and Data Unit
  (HDU) pair from a file and returns a list with data, header, and axis
  information.   Files with more complicated structures or
  isolated header units can be read with an appropriate combination of
  \code{readFITSheader}, \code{readFITSarray}, and
  \code{readFITSbintable}.  See the Example section for
  \code{readFITSimage} for a step-by-step
  example that includes file opening and closing.

  \var{phdu} = 0 forces a read of a secondary HDU in a
  fairly pathological case (specifically: indicating
  a read of a secondary HDU by setting NAXISn = 0 for n > 1, but NAXIS != 0
  and NAXIS1 != 0.)  This does not seem to be forbidden by the FITS
  standard but would be unlikely coding.

  \var{...} passes additional values for file reading.  The only use at
  present is to pass values to \var{fixHdr} in \code{\link{readFITSheader}}.
  \var{fixHdr} attempts to fix headers with non-printing characters,
  either by removing them with \code{fixHdr = 'remove'}, reading further
  into the file until 2880 valid characters are present, or by
  substituting spaces for non-printing characters with \code{fixHdr =
  'substitute'}.  This option should be used with caution, as non-printing
  characters should not be in the header in the first place, so this
  option may or may not corrupt the following data.  The default is
  \code{fixHdr = 'none'}.  Partial matching is allowed.

  Binary table bit, complex, and array
  descriptor data types are not implemented in this release due to a lack
  of examples for testing.
}

\value{
  Return values from \code{readFITS} are in a list. Depending on the
  data type, list entries are:
  \item{header}{Vector with input file header up to END statement.}
  \item{hdr}{Vector with \code{keyword = value} pairs parsed from header.}
  \item{imDat}{Data array (image).}
  \item{axDat}{Data frame with axis scaling and labels (image).}
  \item{col}{Data from each column (bintable).}
  \item{colNames}{Vector of column names, TTYPEn FITS variable (bintable).}
  \item{colUnits}{Vector of column units, TUNITn FITS variable (bintable).}
  \item{TNULLn}{Vector of undefined value definitions, FITS variable
        (bintable).}
  \item{TSCALn}{Vector of multipliers for scaling, FITS variable (bintable).}
  \item{TZEROn}{Vector of zeros for scaling, FITS variable (bintable).}
  \item{TDISPn}{Vector of format information, FITS variable (bintable).}

}

\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{http://fits.gsfc.nasa.gov/}}

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
\emph{fv} and \emph{SAOImage DS9} are in active development
with support for
unix, Windows, and Mac OS-X operating systems, and are available at no cost.

See \code{readFrameFromFITS} to read a binary table
directly into an \R data frame.

}

\seealso{
\code{\link{readFITSarray}},
\code{\link{readFITSbintable}},
\code{\link{readFITSheader}},
\code{\link{readFrameFromFITS}},
\code{\link{modifyHeader}},
\code{\link{image}}, \code{\link{par}}
}
\examples{
require(FITSio)

### Image example
## Make test image with axis information, write to disk
Z <- matrix(1:15, ncol = 3)
writeFITSim(Z, file = "test.fits", c1 = "Test FITS file",
            crpix = c(1,1), crvaln = c(10, 100), cdeltn = c(8, 2),
            ctypen = c("Distance", "Time"),
            cunitn = c("Furlongs", "Fortnights"))
## Read back in and display
X <-  readFITS("test.fits")
ax1 <- axVec(1, X$axDat)          # Make axis vector for image
ax2 <- axVec(2, X$axDat)
xlab <- X$axDat$ctype[1]
ylab <- paste(X$axDat$ctype[2], " [", X$axDat$cunit[2], "]", sep = "")
image(ax1, ax2, X$imDat, xlab = xlab, ylab = ylab)
str(X)
X$axDat                           # Display data frame with axis data
X$hdr[1:10]                       # Header sample
X$hdr[which(X$hdr=="BITPIX")+1]   # BITPIX value from header
unlink("test.fits")
## No axis scale markings
image(X$imDat, xlab = xlab, ylab = ylab, xaxt = "n", yaxt = "n")

### Binary table examples
## Bintable with one row and differently multi-dimensional columns
## Either download example file from
## <http://fits.gsfc.nasa.gov/fits_samples.html>
## and use
\dontrun{filename <- "IUElwp25637mxlo.fits"}
## or, for local example use
filename <- system.file("fitsExamples", "IUElwp25637mxlo.fits",
                        package = "FITSio")

Y <- readFITS(filename)
## Look at contents
str(Y)
Y$colNames
str(Y$col)
Y$hdr[which(Y$hdr=="BITPIX")+1]   # BITPIX value from header
plot(Y$col[[5]], ylab = "Value", main = Y$colNames[5], type = "l")

### Simple flat file example
filename <- system.file("fitsExamples", "2008_03_11_203150.fits",
                        package = "FITSio")
Z <- readFITS(filename)
str(Z)
Z$colNames
str(Z$col)
attach(Z)
xc <- which(colNames == "DMJD")
yc <- which(colNames == "TiltX")
xlab <- paste(colNames[xc], " [", colUnits[xc], "]", sep = "")
ylab <- paste(colNames[yc], " [", colUnits[yc], "]", sep = "")
plot(col[[xc]], col[[yc]], xlab = xlab, ylab = ylab, type = "l")
detach(Z)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
