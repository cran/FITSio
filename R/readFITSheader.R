`readFITSheader` <-
function (zz, maxLines = 5000) 
{
### Function gets FITS header
###
### Takes:
  ## File handle: zz
  ## Maximum number of header lines to read: maxLines
### Returns:
  ## Vector with parsed header
### Requires/Used by:
  ## Requires .fitsHdrParse.r 
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001) 
###
### A. Harris, Univ. MD Astronomy, 3/21/08
###
    foundEnd <- FALSE
    hdr <- character()
    i <- 0
    maxHdrs <- maxLines/36  # max number of header units to read
    while (!foundEnd) {
        ## Each header is a set of 36 80-column card images
        tmp <- .fitsHdrParse(readChar(zz, 2880))
        hdr <- c(hdr, tmp$hdrInfo)
        foundEnd <- tmp$foundEnd
        if (i > maxHdrs) 
            stop("Haven't found END in header after ", maxLines, 
                " header lines")
        i <- i + 1
    }
    return(hdr)
}

`.fitsHdrParse` <-
function (hdr_dat) 
{
### Function parses FITS header
### (preceeded by . to hide it from users in package build)
###   
### Takes:
  ## hdr_dat variable from readFITSheader.r
### Returns:
  ## Character vector: hdrInfo
  ## Logical: foundEnd indicating end of header
### Requires/Used by:
  ## Used by readFITSheader.r
###
### A. Harris, Univ. MD Astronomy, 3/20/08
###
    ## Break header data into card images,
    ## then separate comments and break at keyword = value
    num <- 36      # number of card images
    z <- character(num)
    for (i in 1:num) {
        start <- (i - 1) * 80 + 1
        end <- start + 79
        z[i] <- strsplit(strsplit(substr(hdr_dat, start, end), 
            "/")[[1]][1], "=")
    }
    ## Convert to simple vector
    z <- unlist(z)
    ## Strip spaces and single quotes
    for (i in 1:length(z)) {
        z[i] <- gsub(" ", "", z[i])
        z[i] <- gsub("'", "", z[i])
    }
    ## Check for end in header
    foundEnd <- as.logical(length(which(z == "END")))
    ## Return 
    list(hdrInfo = z, foundEnd = foundEnd)
}

