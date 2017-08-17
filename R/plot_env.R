

#' Quick changes to the plotting environment
#'
#' @param a number of rows (passed to mfrow)
#' @param b number of columns (passed to mfrow)
#' @param brewer_n number of colours to retrieve from RColorBrewer palette
#' @param
mypar <- function(a = 1, b = 1, brewer_n = 9, brewer_name = "Set1",
                  cex.lab = 1, cex.main = 1.2, cex.axis = 1,
                  mar = c(2.5, 2.5, 1.6, 1.1), mgp = c(1.5, 0.5, 0), ...) {

    par(mar = mar, mgp = mgp, cex.lab = cex.lab, cex.main = cex.main,
        cex.axis = cex.axis)
    par(mfrow = c(a, b), ...)
    palette(RColorBrewer::brewer.pal(brewer_n, brewer_name))

}

## next one
