readFITSbintable <-
function (zz, hdr)
{
### Reader for FITS BINTABLE arrays
###
### Takes:
  ## File handle: zz
  ## Parsed header vector: hdr
### Returns:
  ## Data for each column: col
  ## Parsed header vector: hdr
  ## Column names: colNames
  ## Column units: colUnits
  ## TNULLn, TSCALn, TZEROn, TDISPn as defined by FITS standard
### Requires/Used by:
  ## Requires readFITSheader.r
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001)
###
### A. Harris, Univ. MD Astronomy, 4/22/08
### Eliminate blanks in TTFORM etc. variables for new header parser
### AH 7/10/09
### Fixed apparent problem in padding calculation and read, 9/28/10 AH
### Added "K" TFORMn; reads and records high and low 32-bit words separately,
###   returning vector with double length.  Package int64 may be useful to
###   reconstruct 64-bit word.  10/26/2012 AH
###
### Known incomplete:
###  * Data types X, C, M, P unimplemented (bit, complex, double complex,
###      array descriptor)
###
### Updated for new header handling, 12/30/12 AH
###
### Added reads for 16X and 32X bit arrays, updated error messaging.
###   2020.11.14 AH
###
### Added reads for 8X and 24X bit arrays, updated error messaging.
###   2020.11.21 AH
###
###


    ## Parse header if full header is supplied instead of parsed version
    if (nchar(hdr[1])==80) hdr <- parseHdr(hdr)

    ## Determine array dimensions
    naxis1 <- as.numeric(hdr[which(hdr == "NAXIS1") + 1])   # bytes per row
    naxis2 <- as.numeric(hdr[which(hdr == "NAXIS2") + 1])   # rows
    tfields <- as.numeric(hdr[which(hdr == "TFIELDS") + 1]) # columns
    ## Get additional memory allocation information
    tmp <- hdr[which(hdr == "PCOUNT") + 1]
    pcount <- ifelse(length(tmp) != 1, 0, as.numeric(tmp))
    if (pcount != 0)
        warning("pcount must be 0 in a bintable")
    tmp <- hdr[which(hdr == "GCOUNT") + 1]
    gcount <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
    if (gcount != 1)
        warning("gcount must be 1 in a bintable")
    ## Put array information into vectors
    TFORMn <- character(tfields)
    TTYPEn <- character(tfields)  # field names
    TUNITn <- character(tfields)
    TNULLn <- integer(tfields)
    TSCALn <- numeric(tfields)
    TZEROn <- numeric(tfields)
    TDISPn <- character(tfields)
    THEAPn <- integer(tfields)
    TDIMn <- character(tfields)
    for (i in 1:tfields) {
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TFORM", i, sep = "")) +
            1])
        TFORMn[i] <- tmp  # strip out spaces
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TTYPE", i, sep = "")) +
            1])
        TTYPEn[i] <- ifelse(length(tmp) != 1, "", tmp)
        tmp <-gsub(" ", "",  hdr[which(hdr == paste("TUNIT", i, sep = "")) +
            1])
        TUNITn[i] <- ifelse(length(tmp) != 1, "", tmp)
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TNULL", i, sep = "")) +
            1])
        TNULLn[i] <- ifelse(length(tmp) != 1, NA, as.numeric(tmp))
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TSCAL", i, sep = "")) +
            1])
        TSCALn[i] <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TZERO", i, sep = "")) +
            1])
        TZEROn[i] <- ifelse(length(tmp) != 1, 0, as.numeric(tmp))
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TDISP", i, sep = "")) +
            1])
        TDISPn[i] <- ifelse(length(tmp) != 1, "", tmp)
        tmp <- gsub(" ", "", hdr[which(hdr == paste("THEAP", i, sep = "")) +
            1])
        THEAPn[i] <- ifelse(length(tmp) != 1, NA, as.numeric(tmp))
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TDIM", i, sep = "")) +
            1])
        TDIMn[i] <- ifelse(length(tmp) != 1, "", tmp)
    }

    ## Work out formats and set up storage for each column.
    bsize <- integer(tfields)
    btype <- integer(tfields)
    bsign <- logical(tfields)
    mult <- integer(tfields)
    col <- vector("list", tfields)
    for (i in 1:tfields) {
        ## Strip spaces, then separate multiplier and type.
        nc <- nchar(TFORMn[i])
        tmp <- substr(TFORMn[i], 1, nc - 1)
        mult[i] <- ifelse(tmp == "", 1, as.numeric(tmp))
        ## Format; btype: 1 character, 2 logical, 3 integer,
        ## 4 numeric, 5 complex, 6 64-bit integer, 7 24-bit bit array
        form <- tolower(substr(TFORMn[i], nc, nc))
        switch(form, l = {  # logical
            bsize[i] <- 1
            btype[i] <- 2
            bsign[i] <- FALSE
        }, b = {            # 8-bit unsigned int
            bsize[i] <- 1
            btype[i] <- 3
            bsign[i] <- FALSE
        }, i = {            # 16-bit signed int (single)
            bsize[i] <- 2
            btype[i] <- 3
            bsign[i] <- TRUE
        }, j = {            # 32-bit signed int (double)
            bsize[i] <- 4
            btype[i] <- 3
            bsign[i] <- TRUE
        }, k = {            # 64-bit signed int (quadruple)
            bsize[i] <- 8
            btype[i] <- 6
            bsign[i] <- TRUE
        }, a = {            # 8-bit character
            bsize[i] <- 1
            btype[i] <- 1
        }, e = {            # 32-bit float (single)
            bsize[i] <- 4
            btype[i] <- 4
        }, d = {            # 64-bit float (double)
            bsize[i] <- 8
            btype[i] <- 4
            bsign[i] <- FALSE
        }, x = {            # some bit arrays read as integers
            if (mult[i] == 8) {
                btype[i] <- 3
                bsize[i] <- 1
                mult[i] <- 1 # 1 byte
            } else if (mult[i] == 16) {
                btype[i] <- 3
                bsize[i] <- 2
                mult[i] <- 1
            } else if (mult[i] == 24) {
                btype[i] <- 7
                mult[i] <- 3 # 3 bytes
            } else if (mult[i] == 32) {
                btype[i] <- 3
                bsize[i] <- 3
                mult[i] <- 1
            } else stop("Unsupported length in bit array \n")
            bsign[i] <- FALSE  # unsigned
        }, stop('Incompatible TFORMn **', toupper(form),
                '** (C, M, P, arbitrary X lengths not yet implemented)\n'))

        ## Set up storage arrays: rows = number of table rows, columns =
        ## multiplier (depth) of the cells in each column.  Characters are an
        ## exception since they return strings.  64-bit integers need two
        ## 32-bit words per value.
        col[[i]] <- switch(btype[i],
                           #1, char
                           array("", dim = c(naxis2, 1)),
                           #2, logical
                           array(FALSE, dim = c(naxis2, mult[i])),
                           #3, integer
                           array(NA, dim = c(naxis2, mult[i])),
                           #4, numeric
                           array(NA, dim = c(naxis2, mult[i])),
                           #5, complex
                           array(NA, dim = c(naxis2, mult[i])),
                           #6, 64-bit integer, 2 32-bit integers
                           array(NA, dim = c(2*naxis2, mult[i])),
                           #7, 24-bit bit array, 32-bit integer
                           array(NA, dim = c(naxis2, 1))
                           )
    }

    ## Read data, row by row
    for (i in 1:naxis2) {
        for (j in 1:tfields) {
            if (btype[j] <= 5) {
                col[[j]][i, ] <- switch(btype[j],
                                 # 1, char
                                 readChar(zz, nchars = mult[j]),
                                 # 2, logical
                                 readChar(zz, nchars = mult[j]),
                                 # 3, integer
                                 readBin(zz, what = integer(), n = mult[j],
                                         size = bsize[j], signed = bsign[j],
                                         endian = "big"),
                                 # 4, numeric
                                 readBin(zz, what = numeric(), n = mult[j],
                                         size = bsize[j],
                                         endian = "big"),
                                 # 5, complex
                                 readBin(zz, what = complex(), n = mult[j],
                                         size = bsize[j],
                                         endian = "big")
                   )
            }
            if (btype[j] == 6) {
                # 64-bit signed integer, 2 32-bit integers
                col[[j]][i, ] <- readBin(zz, what = integer(),
                                    n = mult[j], size = 4, signed = bsign[j],
                                    endian = "big")
                 col[[j]][i+naxis2, ] <- readBin(zz, what = integer(),
                                    n = mult[j], size = 4, signed = bsign[j],
                                    endian = "big")
            }
            if (btype[j] == 7) {
                # 8 and 24-bit bit arrays, 32-bit integer
                if (bsize[j]==1) { # single byte
                    col[[j]][i, 1] <- as.integer(readBin(zz, "raw",n = 1,
                                                 size = 1, endian = "big"))
                    } else { # 24-bit
                        tmp <- as.integer(readBin(zz, "raw", n = 3, size = 1,
                                                  endian = "big"))
                        col[[j]][i, 1] <- tmp[3]
                        col[[j]][i, 1] <- bitwOr(col[[j]][i, 1],
                                                 bitwShiftL(tmp[2], 8))
                        col[[j]][i, 1] <- bitwOr(col[[j]][i, 1],
                                                 bitwShiftL(tmp[1], 16))
                    }
            }
        }
    }

    ## Finish reading block
    nbyte <- naxis1 * naxis2
    nbyte <- ifelse(nbyte%%2880 == 0, 0, 2880 - nbyte%%2880)
    tmp <- readBin(zz, 'raw', nbyte)

    ## Clean up before returning
    for (i in 1:tfields) {
        ## Apply scaling and offset where appropriate
        if (btype[i] >= 3 && (TSCALn[i] != 1
            || TZEROn[i] != 0)) {
            col[[i]] <- (col[[i]]) * TSCALn[i] + TZEROn[i]
        }
        ## Convert 1D arrays to vectors for easier plotting
        if (nrow(col[[i]]) == 1 || ncol(col[[i]]) == 1
            || btype[i] <= 2) {
            col[[i]] <- as.vector(col[[i]])
        }
    }

    ## Return data list: data for each column plus ancillary data.
    list(col = col, hdr = hdr, colNames = TTYPEn, colUnits = TUNITn,
        TNULLn = TNULLn, TSCALn = TSCALn, TZEROn = TZEROn, TDISPn = TDISPn)
}

