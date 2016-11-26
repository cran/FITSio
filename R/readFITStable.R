readFITStable <-
function (zz, hdr)
{
### Reader for FITS TABLE (ASCII) arrays
###
### Takes:
  ## File handle: zz
  ## Parsed header vector: hdr
### Returns:
  ## Data for each column: col
  ## Parsed header vector: hdr
  ## Column names: colNames
  ## Column units: colUnits
  ## TNULLn, TSCALn, TZEROn: as defined by FITS standard
### Requires/Used by:
  ## Requires readFITSheader.r
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001)
###
### A. Harris, Univ. MD Astronomy, 11/24/2016

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
        warning("pcount must be 0 in a TABLE")
    tmp <- hdr[which(hdr == "GCOUNT") + 1]
    gcount <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
    if (gcount != 1)
        warning("gcount must be 1 in a TABLE")
    ## Put array information into vectors
    TBCOLn <- numeric(tfields)
    TFORMn <- character(tfields)
    TSCALn <- numeric(tfields)
    TZEROn <- numeric(tfields)
    TNULLn <- numeric(tfields)
    TTYPEn <- character(tfields)  # field names
    TUNITn <- character(tfields)
    for (i in 1:tfields) {
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TBCOL", i, sep = "")) +
            1])
        TBCOLn[i] <-  as.numeric(tmp)  
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TFORM", i, sep = "")) +
            1])
        TFORMn[i] <- tmp  
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TSCAL", i, sep = "")) +
            1])
        TSCALn[i] <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TZERO", i, sep = "")) +
            1])
        TZEROn[i] <- ifelse(length(tmp) != 1, 0, as.numeric(tmp))
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TNULL", i, sep = "")) +
            1])
        TNULLn[i] <- ifelse(length(tmp) != 1, NA, as.numeric(tmp))
        tmp <- gsub(" ", "", hdr[which(hdr == paste("TTYPE", i, sep = "")) +
            1])
        TTYPEn[i] <- ifelse(length(tmp) != 1, "", tmp)
        tmp <-gsub(" ", "",  hdr[which(hdr == paste("TUNIT", i, sep = "")) +
            1])
        TUNITn[i] <- ifelse(length(tmp) != 1, "", tmp)
    }

    ## Read data and parse, probably clearing out extra characters,
    ## put result into data frame DF

    DF <- readBin(zz, 'raw', naxis1*naxis2)
    DF <- rawToChar(DF)
    DF <- gsub('[\n]', ' ', DF)  # make some control characters FORTRAN-readable
    tt <- file('Rtmp.txt', open='w')
    writeLines(DF, tt)
    close(tt)

    ## work out widths
    widths <- floor(as.numeric(sapply(strsplit(TFORMn, '[aifedgAIFEDG]'), '[[', 2)))
    ## work out columns to skip
    skips <- character(tfields)
    skips[1] <- paste('X', TBCOLn[1]-1, sep='')  # first column
    for (i in 2:tfields) {
        tmp <- TBCOLn[i] - (TBCOLn[i-1] + widths[i-1])
        skips[i] <- paste('X', tmp, sep='')
    }
    ## assemble format vector, interleave skips and reads
    form <- character(2*tfields-1)
    for (i in 1:tfields) {
        form[2*i-1] <- skips[i]
        form[2*i] <- TFORMn[i]
    }
    ## eliminate zero skips
    if (any(form=='X0')) form <- form[-which(form=='X0')]
    ## formatted read of disk file 
    DF <- read.fortran('Rtmp.txt', format=form, strings=FALSE)
    names(DF) <- TTYPEn
    unlink('Rtmp.txt')
    
    ## Finish reading block
    nbyte <- naxis1 * naxis2
    nbyte <- ifelse(nbyte%%2880 == 0, 0, 2880 - nbyte%%2880)
    tmp <- readBin(zz, 'raw', nbyte)

    ## Return data list: data for each column plus ancillary data.
    list(hdr = hdr, DF=DF, colNames = TTYPEn, colUnits = TUNITn,
         TNULLn = TNULLn, TSCALn = TSCALn, TZEROn = TZEROn)
}

