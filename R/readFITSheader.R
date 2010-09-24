`readFITSheader` <-
function (zz, maxLines = 5000, fixHdr = 'none')
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
  ## Added header length check 9/22/10 AH
  ## Added header fixes (remove/replace) 9/25/2010 AH
###
    foundEnd <- FALSE
    hdr <- character()
    i <- 0
    maxHdrs <- maxLines/36  # max number of header units to read
    while (!foundEnd) {
        ## Each header is a set of 36 80-column card images
        switch(pmatch(fixHdr, c('none', 'remove', 'substitute'), nomatch=4),
               { # header is ok
                    inpString <- readChar(zz, 2880)
                    if (nchar(inpString) != 2880) {
                        txt <- paste('*** Header problem:', nchar(inpString),
                                     'characters instead of 2880; try option fixHdr *** \n')
                        close(zz)
                        stop(txt)
                    }
                }, { # substitute blanks for non-printing characters from header
                   nbytes <- 2880
                   inpBin <- raw()
                   while (nbytes != 0) {
                       inpBin <- c(inpBin, readBin(zz, 'raw', nbytes))
                       idxBad <- which(inpBin <= 0x1F | inpBin == 0x7F)
                       nbytes <- length(idxBad)
                       if (nbytes > 0) {
                           inpBin <- inpBin[-idxBad]
                           txt <- paste('*** Removed', length(idxBad),
                                        'non-printing characters in header ***\n')
                           cat(txt)
                       }
                       inpString <- rawToChar(inpBin)
                   }
                }, { # substitute spaces for non-printing characters
                    inpBin <- readBin(zz, 'raw', 2880)
                    idxBad <- which(inpBin <= 0x1F | inpBin == 0x7F)
                    nbytes <- length(idxBad)
                    if (nbytes > 0) {
                        inpBin[idxBad] <- as.raw(0x20)  # substitute space
                        txt <- paste('*** Substituted', nbytes,
                                     'spaces for non-printing characters ***\n')
                        cat(txt)
                    }
                    inpString <- rawToChar(inpBin)
                }, {
                    txt <- '*** fixHdr must be one of: none, remove, or substitute ***\n'
                    close(zz)
                    stop(txt)
                }

        )
        tmp <- .fitsHdrParse(inpString)
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
