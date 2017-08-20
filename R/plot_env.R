#' Quick changes to the plotting environment
#'
#' @param a number of rows (passed to mfrow)
#' @param b number of columns (passed to mfrow)
#' @param brewer_n number of colours to retrieve from RColorBrewer palette
#' @param brewer_name RColorBrewer palette
#' @param ... further arguments passed to par
#' @return By default only a new colour palette
mypar <- function(a = 1, b = 1, brewer_n = 9, brewer_name = "Set1", ...) {

    par(mfrow = c(a, b), ...)
    palette(RColorBrewer::brewer.pal(brewer_n, brewer_name))

}

