\name{readFITSheader}
\alias{readFITSheader}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a FITS header}
\description{
  Read a FITS header from an open connection to a FITS file.
}
\usage{
readFITSheader(zz, maxLines = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zz}{file handle; see Example.}
  \item{maxLines}{maximum number of header lines to read; see Details.}
}
\details{
  \code{readFITSheader} reads the data from the header part of a FITS
  Header and Data Unit.  In addition to general header information, it
  provides parameters needed to read image and binary table files by
  functions \code{readFITSarray} and \code{readFITSbintable}.  A header
  unit may exist without a corresponding data unit, for instance to
  carry additional coordinate information.

  The \var{maxLines} variable limits the number of header lines
  \code{readFITSheader} will read to prevent endless reading if the
  header is flawed and the END
  statement is missing.  The function generates a message an halts when the
  number of reads exceeds  \var{maxLines}.  Increase the value as
  needed for very large headers.

  The header vector does not have an easy format for
  people to read (graphical FITS viewers like \emph{fv} and
  \emph{SAOImage DS9} are good for this),
  but is designed for further processing
  by \R.  The header vector has a kind of
  \emph{keyword value keyword value \ldots}
  format, where keywords without values are simply followed by the next
  keyword.  This means a seach for a keyword will give the corresponding
  value as the next element in the vector; see the Examples.    
}
\value{
  \var{hdr} is a character vector \var{hdr}
  containing the header information in a format that is easy for
  \R to parse further.  See Details.
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

\seealso{
  \code{\link{readFITS}}, \code{\link{readFITSarray}},
  \code{\link{readFITSbintable}}, \code{\link{file}}
}

\examples{
require("FITSio")
## Make test image with axis information, write to disk
Z <- matrix(1:15, ncol = 3)
writeFITSim(Z, file = "test.fits", c1 = "Test FITS file",
            crpix = c(1,1), crvaln = c(10, 100), cdeltn = c(8, 2),
            ctypen = c("Distance", "Time"),
            cunitn = c("Furlongs", "Fortnights"))
## Read back in 
## Open file, read header and array, close file and delete.
    zz <- file(description = "test.fits", open = "rb")
    hdr <- readFITSheader(zz)
    D <- readFITSarray(zz, hdr) 
    close(zz)
    hdr[1:10]                       # Header sample
    hdr[which(hdr=="BITPIX")+1]   # BITPIX value from header
    unlink("test.fits")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}

