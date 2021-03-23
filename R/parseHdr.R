parseHdr <-
function(headerName) {
## Take header card images and produce a vector with elements keyword, value
##
## A. Harris 2012.10.30, 2013.01.28, 2016.12.18

    # identify keyword=value pairs by checking for '=' as 9th 
    # character (strict FITS encoding), keep only those pairs
    headerName <- headerName[which(substr(headerName, 9, 9)=='=')]

    # make hdr array
    len <- length(headerName)
    hdr <- character(2*len)
    
    # separate keyword, value/notes
    for (i in 1:len) {
        hdr[2*i-1] <- substr(headerName[i], 1, 8)
        hdr[2*i] <- substr(headerName[i], 10, 80)
    }
    
    # find and extract string keywords
    idx <- grep("^ *'", hdr)
    for (i in idx) {
        # replace doubled single quotes with improbable dummy string
        hdr[i] <- gsub("''", 'aAlJ2fZ47xx', hdr[i])
        # extract keyword: string between single quotes
        hdr[i] <- strsplit(hdr[i], "'")[[1]][2]
        # replace improbable dummy string with doubled single quotes
        hdr[i] <- gsub('aAlJ2fZ47xx', "''", hdr[i])
    }
    
    # find numerical keywords, including negative, get rid of notes
    idx <- grep("^ *-*[[:digit:]]", hdr)
    for (i in idx) {
        hdr[i] <- strsplit(hdr[i], '/')[[1]][1]
    }

    # find logical keywords, get rid of notes
    idx <- grep("^ *[TFtf]", hdr)
    for (i in idx) {
        hdr[i] <- strsplit(hdr[i], '/')[[1]][1]
    }
        
    # eliminate leading and trailing spaces
    for (i in 1:length(hdr)) {
        hdr[i] <- sub('^ *', '', hdr[i])
        hdr[i] <- sub(' *$', '', hdr[i])
    }

    # return parsed header
    hdr
}
