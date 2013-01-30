\name{parseHdr}
\alias{parseHdr}
\title{
Parse FITS header
}
\description{
Parse FITS header to extract \code{keyword = variable} pairs.  Mostly
useful for finding values given keyword name.
}
\usage{
parseHdr(headerName)
}
\arguments{
  \item{headerName}{
Header card images (vector)
}
}
\details{
Uses strict FITS convention to identify  \code{keyword = variable}
pairs: a string '= ' in card image columns 9 and 10.
}
\value{
Vector with elements \code{keyword, variable} for all pairs.
}

\author{
A. Harris
}
\note{
Function eliminates leading and trailing spaces.
}

\seealso{
\code{\link{readFITSheader}}, \code{\link{modifyHeader}}
}
\examples{
header <- newKwv('KEYWORD', 'VALUE', 'NOTE')
header <- addKwv('test1', 'plot size', header=header)
header <- addKwv('test2', 4294.95397809807, 'number', header=header)
header <- addKwv('test3', 4.29495397809807e50, 'big number', header=header)
header <- addKwv('test4', -4.29495397809807e50, 'big number', header=header)
parseHdr(header)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{files}
