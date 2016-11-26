\name{readFITStable}
\alias{readFITStable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a FITS ASCII table}
\description{Read a FITS ASCII table from an open connection to a
FITS file.}

\usage{
readFITStable(zz, hdr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zz}{File handle; see Example.}
  \item{hdr}{Header card images, raw or parsed.}
}
\details{
  \code{readFITStable} reads the data from the ASCII table part of a FITS
  Header and Data Unit (hdu) containing ASCII table data.  The header must be
  read first by \code{\link{readFITSheader}}; either this header or the
  parsed version from \code{\link{parseHdr}} are valid for the \var{hdr}
  variable.  ASCII tables are multi-column files with one or more rows.  

}
\value{
  \item{hdr}{Vector with parsed header.}
  \item{DF}{Data frame containing table data.}
  \item{colNames}{Vector of column names, TTYPEn FITS variable.}
  \item{colUnits}{Vector of column units, TUNITn FITS variable.}
  \item{TNULLn}{Vector of undefined value definitions, FITS variable.}
  \item{TSCALn}{Vector of multipliers for scaling, FITS variable.}
  \item{TZEROn}{Vector of zeros for scaling, FITS variable.}
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
  \code{\link{readFITS}}, \code{\link{readFITSheader}},
  \code{\link{readFITSarray}}, \code{\link{file}}
}

\examples{
require(FITSio)

filename <- system.file("fitsExamples", "vizier.fits",
                        package = "FITSio")

## Simple read
D <- readFITS(filename)

## Explicit read: open file, read header and table, close file.
zz <- file(description = filename, open = "rb")
header0 <- readFITSheader(zz) # read primary header
header <- readFITSheader(zz) # read extension header
D <- readFITStable(zz, header)
close(zz)

## Either way, look at contents
str(D)
str(D$DF)
str(header)
str(parseHdr(header))
D$hdr[which(D$hdr=="BITPIX")+1]   # BITPIX value from header
D$colNames

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}

