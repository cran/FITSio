\name{readFITSheader}
\alias{readFITSheader}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a FITS header}
\description{
  Read a FITS header from an open connection to a FITS file.
}
\usage{
readFITSheader(zz, maxLines = 5000, fixHdr = 'none')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zz}{file handle; see Example.}
  \item{maxLines}{maximum number of header lines to read; see Details.}
  \item{fixHdr}{deal with non-printing characters in header; see Details.}
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

  \var{fixHdr} attempts to fix headers with non-printing characters,
  either by removing them with \code{fixHdr = 'remove'}, reading further
  into the file until 2880 valid characters are present, or by
  substituting spaces for non-printing characters with \code{fixHdr =
  'substitute'}.  This option should be used with caution, as non-printing
  characters should not be in the header in the first place, so this
  option may or may not corrupt the following data.  The default is
  \code{fixHdr = 'none'}.  Partial matching is allowed.

  The header vector does not have an easy format for
  people to read (graphical FITS viewers like \emph{fv} and
  \emph{SAOImage DS9} are good for this),
  but is designed for further processing
  by \R.  The header vector has a kind of
  \emph{keyword value keyword value \ldots}
  format, where keywords without values are simply followed by the next
  keyword.  This means a search for a keyword will give the corresponding
  value as the next element in the vector; see the Examples.
}
\value{
  \var{hdr} is a character vector \var{hdr}
  containing the header information in a format that is easy for
  \R to parse further.  See Details.
}
\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{https://fits.gsfc.nasa.gov/}
}

\author{Andrew Harris}

\note{
  Graphical FITS viewers such as \emph{fv}
  (\url{https://heasarc.gsfc.nasa.gov/ftools/fv/}) and \emph{SAOImage DS9}
(\url{http://ds9.si.edu/})
have excellent facilities for displaying FITS data, headers, and file
structure.  Having one or more graphical viewers available will
prove extremely
useful
for working with FITS
files, even when the data are read into \R for further processing.
\emph{fv} and \emph{SAOImage DS9} are in active development with support for
unix, Windows, and Mac OS-X operating systems, and are available at no cost.
}

\seealso{
  \code{\link{parseHdr}}, the usual complement to \code{readFITSheader};
  \code{\link{readFITS}}, \code{\link{readFITSarray}},
  \code{\link{readFITSbintable}}, \code{\link{file}}
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
## Read back in
## Open file, read header and array.
    zz <- file(description = filename, open = "rb")
    header <- readFITSheader(zz)
    hdr <- parseHdr(header)
    D <- readFITSarray(zz, hdr)
    close(zz)
    hdr[1:10]                       # Header sample
    hdr[which(hdr=="BITPIX")+1]   # BITPIX value from header

## Clean up files to avoid clutter
unlink(filename)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}

