addHistory <-
function(history, headerName) {
# Add comment line card image to FITS header
#
# A. Harris 2012.10.13
    # check for valid inputs, stop otherwise
    if (!is.character(history)) stop('*** History must be a string ***')
    history <- strtrim(history, 70)

    # modify header
    headerName <- c(headerName, sprintf('HISTORY   %-70s', history))

    # return modified structure
    headerName
}
