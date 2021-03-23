\name{makeFITSimHdr}
\alias{makeFITSimHdr}
\title{
Generate header for FITS image
}
\description{
Function generates the header for a FITS image.  It determines much of
the header information from the file itself and can merge header
information with error checking and removal of conflicting reserved keywords.
}
\usage{
makeFITSimHdr(X, primaryhdu = TRUE, type = "double", c1 = NA, c2 = NA,
bscale = 1, bzero = 0, crpixn = NA, crvaln = NA, cdeltn = NA, ctypen =
NA, cunitn = NA, axDat = NA, header = "")
}
\arguments{
  \item{X}{Multi-dimensional numeric data array; see Details.}
  \item{primaryhdu}{Logical: TRUE for stand-alone file, FALSE if
    subordinate header.}
  \item{type}{Type to write: single or double precision.}
  \item{c1}{Character string comment line for header.}
  \item{c2}{Character string comment line for header.}
  \item{bscale}{Global scaling factor, FITS standard meaning.}
  \item{bzero}{Global shift, FITS standard meaning.}
  \item{crpixn}{Vector of reference pixel numbers for axes, FITS
    standard meaning.}
  \item{crvaln}{Vector of values at reference pixels, FITS standard meaning.}
  \item{cdeltn}{Vector of axis increments per pixel, FITS standard meaning.}
  \item{ctypen}{String vector of descriptive labels for axis, FITS
    standard meaning.}
  \item{cunitn}{String vector of physical units for axis, FITS standard
    meaning.}
  \item{axDat}{Data frame with axis data, see details.}
  \item{header}{Optional header of 80-character card images.}
}

\value{
String with format and length matching FITS header standard for an image.
}
\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{https://fits.gsfc.nasa.gov/}
}
\author{
Andrew Harris
}

\seealso{
\code{\link{writeFITSim}}
}
\examples{
X <- matrix(1:15, ncol = 3)
# extra header lines (optional)
header <- newKwv('KEYWORD', 'VALUE', 'NOTE') # initialize header
header <- addComment('Add these lines to auto-generated header', header=header)
header <- delKwv('KEYWORD', header=header)  # kill first line
header <- addKwv('test1', 'plot size', header=header)
header <- addKwv('test2', 4294.95397809807, 'number', header=header)
header <- addKwv('test3', 4.29495397809807e50, 'big number', header=header)
# make main header
header <- makeFITSimHdr(X, crpixn = c(1,1), crvaln = c(10, 100),
                        cdeltn = c(8, 2), ctypen = c("Distance", "Time"),
                        cunitn = c("Furlongs", "Fortnights"),
                        header = header)
# finish and close out header
tmp <- closeHdr(header) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
