`writeFITSim16i` <-
function (X, file = 'R.fits', ...) 
{
### Utility writes FITS image after converting to single integer words
### for more compact data files
###
### Takes:
  ## Multi-dimensional array: X
  ## Output FITS file: file
  ## Other options for writeFITSim.r (but not bscale or bzero!)
### Returns:
  ## Writes FITS file to disk
### Requires/Used by:
  ## Requires writeFITSim.r
###  
### A. Harris, Univ. MD Astronomy, 4/22/08, 7/24/08
###
    ## Set scale and offset to use maximum range
    scale <- (2^16 - 2)/(max(X) - min(X))
    offset <- (-2^15 + 1)/scale - min(X)
    ## Scale and offset data, convert to integer array
    X <- (X + offset) * scale
    X <- array(as.integer(X), dim = dim(X))
    ## Write data
    writeFITSim(X, file, type = "single", bscale = 1/scale, 
        bzero = -offset, ...)
}

