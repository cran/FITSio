\name{readFITSbintable}
\alias{readFITSbintable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a FITS binary table}
\description{Read a FITS binary table from an open connection to a
FITS file.}

\usage{
readFITSbintable(zz, hdr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zz}{file handle; see Example.}
  \item{hdr}{vector of header information from \code{readFITSheader}.}
}
\details{
  \code{readFITSbintable} reads the data from the binary table part of a FITS
  Header and Data Unit containing binary table data.  The header must be
  read first by \code{readFITSheader} to provide the variable \var{hdr}.
  Binary tables are multi-column files with one or more rows.  Each
  column has an individual data type and number of entries per cell
  (i.e., a cell may contain a scalar or vector).

  Binary table bit, complex, and array
  descriptor data types are not implemented in this release due to a lack
  of examples for testing.
}
\value{
  \item{col}{Data from each column, either a vector or an array.}
  \item{hdr}{Vector with parsed header.}
 \item{colNames}{Vector of column names, TTYPEn FITS variable.}
  \item{colUnits}{Vector of column units, TUNITn FITS variable.}
  \item{TNULLn}{Vector of undefined value definitions, FITS variable.}
  \item{TSCALn}{Vector of multipliers for scaling, FITS variable.}
  \item{TZEROn}{Vector of zeros for scaling, FITS variable.}
  \item{TDISPn}{Vector of format information, FITS variable.}
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
require("FITSio")

## Either download example file from
## <http://fits.gsfc.nasa.gov/fits_samples.html>
## and use
\dontrun{filename <- "IUElwp25637mxlo.fits"}
## or, for local example use
filename <- system.file("fitsExamples", "IUElwp25637mxlo.fits",
                        package = "FITSio")

## Open file, read header and table, close file.
zz <- file(description = filename, open = "rb")
hdr <- readFITSheader(zz) # read primary header
hdr <- readFITSheader(zz) # read extension header
D <- readFITSbintable(zz, hdr)
close(zz)

## Look at contents
str(D)
D$hdr[which(D$hdr=="BITPIX")+1]   # BITPIX value from header
D$colNames
plot(D$col[[5]], ylab = "Value", main = D$colNames[5], type = "l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}

