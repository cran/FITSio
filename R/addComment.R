addComment <-
function(comment, headerName) {
# Add comment line card image to FITS header
#
# A. Harris 2012.10.13
    # check for valid inputs, stop otherwise
    if (!is.character(comment)) stop('*** Comment must be a string ***')
    comment <- strtrim(comment, 70)

    # modify header
    headerName <- c(headerName, sprintf('COMMENT   %-70s', comment))

    # return modified structure
    headerName
}
