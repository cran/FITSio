\name{closeHdr}
\alias{closeHdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format and close FITS header
}
\description{
Function adds END statement and closes FITS header
}
\usage{
closeHdr(headerName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{headerName}{
Header card images (vector)
}
}
\details{
Adds \code{END} to header and pads with spaces to length defined in FITS
standard.
}
\value{
Character vector to write to FITS file as header.
}
\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{http://fits.gsfc.nasa.gov/}
}
\author{
A. Harris
}

\seealso{
\code{\link{modifyHeader}}, \code{\link{makeFITSimHdr}}
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
\keyword{file}

