\name{axVec}
\alias{axVec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate axis vector for image}
\description{
  \code{axVec} generates a vector of axis values from variables
  contained in a FITS image file header.
}
\usage{
axVec(nax = 1, axDat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nax}{index number of the axis.}
  \item{axDat}{data table axis variables produced by
    \code{readFITSarray} (which is also contained in \code{readFITS}).}
}
\details{
  Run once for each axis needed.
}
\value{
  A vector with length equal to the number of pixels along the axis.
  Vector values \eqn{v[i]} for \eqn{i = 1:NAXISn} (\eqn{NAXISn = len},
  the vector length, in
  the axDat data frame) are:
  \deqn{ v[i] = (i - CRPIXn) CDELTn + CRVALn}
  where \eqn{CRPIXn}, \eqn{CDELTn}, and \eqn{CVALn} are the reference
  pixel, pixel increment, and reference pixel value for each axis
  \eqn{n}, following the FITS standard.
}

\references{
Hanisch et al., \emph{Astron.\ Astrophys.} 376, 359-380 (2001)

\url{http://fits.gsfc.nasa.gov/}}

\author{Andrew Harris}

\seealso{\code{\link{readFITS}}, \code{\link{readFITSarray}}}

\examples{
require(FITSio)
Z <- matrix(1:15, ncol = 3)
writeFITSim(Z, file = 'test.fits', c1 = 'Test FITS file',
            crpix = c(1,1), crvaln = c(10, 100), cdeltn = c(8, 2),
            ctypen = c('Distance', 'Time'),
            cunitn = c('Furlongs', 'Fortnights'))
X <-  readFITS('test.fits')
ax1 <- axVec(1, X$axDat)
ax2 <- axVec(2, X$axDat)
xlab <- X$axDat$ctype[1]
ylab <- paste(X$axDat$ctype[2], " [", X$axDat$cunit[2], "]", sep = "")
image(ax1, ax2, X$imDat, xlab = xlab, ylab = ylab)
unlink("test.fits")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
