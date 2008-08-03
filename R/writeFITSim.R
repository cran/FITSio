`writeFITSim` <-
function (X, file = "R.fits", type = "double", bscale = 1, 
    bzero = 0, c1 = NA, c2 = NA, crpixn = NA, crvaln = NA, cdeltn = NA, 
    ctypen = NA, cunitn = NA) 
{
### Function writes FITS image after primary header
###
### Takes:
  ## Multi-dimensional array: X
  ## Output file name: file
  ## Type to write (single or double precision)
  ## Overall scaling and shifting variables: bscale, bzero
  ## Text comment variables: c1 and c2
  ## Axis parameters with usual FITS standard meanings 
### Returns:
  ## Writes FITS file to disk
### Requires/Used by:
  ## Requires .writeFITShdr0.r
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001) 
###
### A. Harris, Univ. MD Astronomy, 4/22/08
###
    ## Number of axes, parse data type (single or double precision)
    naxis <- length(dim(X))
    type <- tolower(substr(type, 1, 1))
    ## Choose integer, float; byte single, double precision based on data
    if (is.integer(X)) {
        switch(type, b = {
            bitpix <- 8
            size <- 1
        }, s = {             
            bitpix <- 16
            size <- 2
        }, d = {
            bitpix <- 32
            size <- 4
        }, stop("Unrecognized data type"))
    }
    else {
        switch(type, s = {
            bitpix <- -32
            size <- 4
        }, d = {
            bitpix <- -64
            size <- 8
        }, stop("Unrecognized data type"))
    }
    ## Open file
    zz <- file(file, "wb")
    ## Write primary header (do not modify bitpix, naxis, or naxisn!)
    hdr0 <- .writeFITShdr0(bitpix = bitpix, naxis = naxis, naxisn = dim(X), 
        bscale = bscale, bzero = bzero, c1 = c1, c2 = c2, crpixn = crpixn, 
        crvaln = crvaln, cdeltn = cdeltn, ctypen = ctypen, cunitn = cunitn)
    writeChar(hdr0, zz, eos = NULL)
    ## Then write data
    writeBin(as.vector(X), zz, size = size, endian = "big")
    ## Then pad rest of record with zeros
    pad <- raw(2880 - (length(as.vector(X)) * size)%%2880)
    writeBin(pad, zz, endian = "big")
    ## Finish up
    close(zz)
}

`.writeFITShdr0` <-
function (bitpix = bitpix, naxis = naxis, naxisn = naxisn, bscale = 1, 
    bzero = 0, c1 = NA, c2 = NA, crpixn = NA, crvaln = NA, cdeltn = NA, 
    ctypen = NA, cunitn = NA) 
{
### Function assembles FITS primary header for images
###    (multi-dimensional arrays)
### (preceeded by . to hide it from users in package build)
###
### Takes:
  ## Most variables have names as defined in FITS reference
  ## Additional comment lines: c1, c2
### Returns:
  ## Header data for writeFITSim.r
### Requires/Used by:
  ## 
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001) 
###
### A. Harris, Univ. MD Astronomy, 4/17/08
###
    ## Make defaults, also if arguments are missing
    if (is.na(crpixn[1])) 
        crpixn[1:naxis] <- 1
    if (is.na(crvaln[1])) 
        crvaln[1:naxis] <- 1
    if (is.na(cdeltn[1])) 
        cdeltn[1:naxis] <- 1
    if (is.na(ctypen[1])) 
        ctypen[1:naxis] <- ""
    if (is.na(cunitn[1])) 
        cunitn[1:naxis] <- ""
    bpad <- sprintf("%80s", " ")
    ## Assemble header string
    txt <- "SIMPLE  =                    T / file conforms to FITS standard"
    hdr <- strtrim(paste(txt, bpad), 80)
    txt <- sprintf("BITPIX  = %20d / number of bits per data pixel", 
        bitpix)
    hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    txt <- sprintf("NAXIS   = %20d / number of data axes", naxis)
    hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    for (i in 1:naxis) {
        lab <- paste("NAXIS", i, sep = "")
        txt <- sprintf("%s  = %20d / length of data axis", lab, 
            naxisn[i])
        hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    }
    txt <- "EXTEND  =                    T / FITS dataset may contain extensions"
    hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    if (!is.na(c1)) {
        txt <- c1
        hdr <- paste(hdr, strtrim(paste("COMMENT  ", txt, bpad), 
            80), sep = "")
    }
    if (!is.na(c2)) {
        txt <- c2
        hdr <- paste(hdr, strtrim(paste("COMMENT  ", txt, bpad), 
            80), sep = "")
    }
    txt <- "COMMENT   Written by writeFITSim.r  ver 1.1"
    hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    txt = "COMMENT   FITS (Flexible Image Transport System) format is defined in"
    hdr = paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    txt <- "COMMENT   Astronomy and Astrophysics, volume 376, page 359 (2001)"
    hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    txt <- sprintf("BSCALE  = %20.9G / overall scaling", bscale)
    hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    txt <- sprintf("BZERO   = %20.9G / overall offset", bzero)
    hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    for (i in 1:naxis) {
        lab <- paste("CRPIX", i, sep = "")
        txt <- sprintf("%s  = %20.9G", lab, crpixn[i])
        hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
        lab <- paste("CRVAL", i, sep = "")
        txt <- sprintf("%s  = %20.9G", lab, crvaln[i])
        hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
        lab <- paste("CDELT", i, sep = "")
        txt <- sprintf("%s  = %20.9G", lab, cdeltn[i])
        hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
        lab <- paste("CTYPE", i, sep = "")
        txt <- sprintf("%s  = '%s'", lab, ctypen[i])
        hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
        lab <- paste("CUNIT", i, sep = "")
        txt <- sprintf("%s  = '%s'", lab, cunitn[i])
        hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    }
    txt <- "END"
    hdr <- paste(hdr, strtrim(paste(txt, bpad), 80), sep = "")
    ## Fill header to 2880 bytes 
    ncards <- nchar(hdr, type = "chars")/80
    ncards <- ifelse(ncards%%36 == 0, 0, (1 - (ncards/36)%%1) * 
        36)
    if (ncards > 0) {
        for (i in 1:ncards) hdr <- paste(hdr, bpad, sep = "")
    }
    ## Done, return header
    return(hdr)
}

