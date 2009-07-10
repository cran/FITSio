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
### Updated to properly parse and preserve values in quotes (but
### deletes leading and trailing blanks in quotes to preserve compatibility
### with existing code)
### AH 7/11/09
###
    ## Break header data into card images,
    ## then separate comments and break at keyword = value
    num <- 36      # number of card images
    z <- character(2*num)  # keyword, value
    i <- 0
    for (j in 1:num) {
        start <- (j - 1) * 80 + 1
        end <- start + 79
        image <- substr(hdr_dat, start, end)
        if (toupper(substr(image, 1, 7)) != 'COMMENT' &&
            toupper(substr(image, 1, 7)) != 'HISTORY' &&
            substr(gsub(" ", "", image), 1, 1) != '/' &&
            gsub(" ", "", image) != "") {

            i <- i+1
            idx <- i*2
            # Parse keyword and value from each line, remove trailing blanks
            # two cases: one with value in quotes (preserve string), one not
            tmp = unlist(strsplit(image, '\''))
            if (length(tmp) == 1) {  # for value not in quotes
                tmp <- unlist(strsplit(strsplit(tmp, '/')[[1]][1], '='))
                z[idx-1] <- gsub(' ', '', tmp[1])
                z[idx] <- gsub(' ', '', tmp[2])
            } else {                # for value in quotes
                tmp2 = unlist(strsplit(tmp[1], '='))[1]
                z[idx-1] <- gsub(' ', '', tmp2)
                z[idx-1] <- gsub("=", "", z[idx-1])
                z[idx] <- sub(' +$', '', tmp[2])
                z[idx] <- sub('^ +', '', z[idx])
            }
        }
    }

    # trim off unused elements at end
    z <- z[1:(i*2)]

    ## Check for end in header
    foundEnd <- as.logical(length(which(z == "END")))
    if (foundEnd) z <- z[1:which(z == "END")]

    ## Return
    list(hdrInfo = z, foundEnd = foundEnd)
}
