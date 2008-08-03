`readFITSbintable` <-
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
###
### Known incomplete:
###  * Data types X, C, M, P unimplemented (bit, complex, double complex, 
###      array descriptor)
###
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
        tmp <- hdr[which(hdr == paste("TFORM", i, sep = "")) + 
            1]
        TFORMn[i] <- tmp
        tmp <- hdr[which(hdr == paste("TTYPE", i, sep = "")) + 
            1]
        TTYPEn[i] <- ifelse(length(tmp) != 1, "", tmp)
        tmp <- hdr[which(hdr == paste("TUNIT", i, sep = "")) + 
            1]
        TUNITn[i] <- ifelse(length(tmp) != 1, "", tmp)
        tmp <- hdr[which(hdr == paste("TNULL", i, sep = "")) + 
            1]
        TNULLn[i] <- ifelse(length(tmp) != 1, NA, as.numeric(tmp))
        tmp <- hdr[which(hdr == paste("TSCAL", i, sep = "")) + 
            1]
        TSCALn[i] <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
        tmp <- hdr[which(hdr == paste("TZERO", i, sep = "")) + 
            1]
        TZEROn[i] <- ifelse(length(tmp) != 1, 0, as.numeric(tmp))
        tmp <- hdr[which(hdr == paste("TDISP", i, sep = "")) + 
            1]
        TDISPn[i] <- ifelse(length(tmp) != 1, "", tmp)
        tmp <- hdr[which(hdr == paste("THEAP", i, sep = "")) + 
            1]
        THEAPn[i] <- ifelse(length(tmp) != 1, NA, as.numeric(tmp))
        tmp <- hdr[which(hdr == paste("TDIM", i, sep = "")) + 
            1]
        TDIMn[i] <- ifelse(length(tmp) != 1, "", tmp)
    }
    ## Work out formats and set up storage for each column.
    bsize <- integer(tfields)
    btype <- integer(tfields)
    bsign <- logical(tfields)
    mult <- integer(tfields)
    col <- vector("list", tfields)
    for (i in 1:tfields) {
        ## Separate multiplier and type.
        nc <- nchar(TFORMn[i])
        tmp <- substr(TFORMn[i], 1, nc - 1)
        mult[i] <- ifelse(tmp == "", 1, as.numeric(tmp))
        ## Format; btype: 1 character, 2 logical, 3 integer,
        ## 4 numeric, 5 complex
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
        }, a = {            # 8-bit character
            bsize[i] <- 1
            btype[i] <- 1
            bsign[i] <- FALSE
        }, e = {            # 32-bit float (single)
            bsize[i] <- 4
            btype[i] <- 4
            bsign[i] <- FALSE
        }, d = {            # 64-bit float (double)
            bsize[i] <- 8
            btype[i] <- 4
            bsign[i] <- FALSE
        }, stop("X, C, M, P not yet implemented \n"))
        ## Set up storage arrays: rows = number of table rows, columns = 
        ## multiplier (depth) of the cells in each column.  Characters are an
        ## exception since they return strings. 
        if (btype[i] == 1) {
            col[[i]] <- array("", dim = c(naxis2, 1))
        }
        else {
            col[[i]] <- array(switch(btype[i], "", FALSE, NA, 
                NA, NA), dim = c(naxis2, mult[i]))
        }
    }
    ## Read data, row by row
    for (i in 1:naxis2) {
        for (j in 1:tfields) {
            if (btype[j] <= 2) {  # character reads
                col[[j]][i, ] <- readChar(zz, nchars = mult[j])
            }
            else {
                what <- switch(btype[j], character(), logical(), 
                  integer(), numeric(), complex())
                col[[j]][i, ] <- readBin(zz, what = what, n = mult[j], 
                  size = bsize[j], signed = bsign[j], endian = "big")
            }
        }
    }
    ## Finish reading block
    nbyte <- naxis1 * naxis2
    nbyte <- ifelse(nbyte%%2880 == 0, 0, (1 - (nbyte/2880)%%1) * 
        2880)
    tmp <- readChar(zz, nbyte)
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
        ## Terminate character strings that end in ascii 0
        if (btype[i] <= 2) {
            lcol <- length(col[[i]])
            txttmp <- character(lcol)
            for (j in 1:lcol) txttmp[j] <-
                strsplit((col[[i]][j]), 0)
            col[[i]] <- unlist(txttmp)
        }
    }
    ## Return data list: data for each column plus ancillary data.  
    list(col = col, hdr = hdr, colNames = TTYPEn, colUnits = TUNITn, 
        TNULLn = TNULLn, TSCALn = TSCALn, TZEROn = TZEROn, TDISPn = TDISPn)
}

