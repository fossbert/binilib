#' Quick changes to the plotting environment
#'
#' @param a number of rows (passed to mfrow)
#' @param b number of columns (passed to mfrow)
#' @param brewer_n number of colours to retrieve from RColorBrewer palette
#' @param brewer_name RColorBrewer palette
#' @param ... further arguments passed to par
#' @return By default only a new colour palette
#' @export
mypar <- function(a = 1, b = 1, brewer_n = 9, brewer_name = "Set1", ...) {

    graphics::par(mfrow = c(a, b), ...)
    grDevices::palette(RColorBrewer::brewer.pal(brewer_n, brewer_name))
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 10))
}

#' Plot quantiles per column for a numeric matrix
#'
#' @param eset numeric matrix typically containing gene expression data
#' @param quants integer of number of top MR to consider
#' @param brewer_name RColorBrewer palette, defaults to Dark2
#' @param ... further arguments passed to matplot
#' @return graph with the respective quantiles across the sample space
#' @export

kaboxplot <- function(eset,
                      quants = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      brewer_name = 'Dark2',
                      ylab = 'Normalized expression',
                      xlab = 'Samples', ...) {

    nq <- length(quants)
    qs <- t(apply(eset, 2, quantile, probs = quants, na.rm = TRUE))
    info <- colMeans(qs)

    oldpar <- par()$mar
    par(mar = rep(4.1, 4), las = 1)
    matplot(qs,
            type = "l",
            lty = 1,
            col = RColorBrewer::brewer.pal(nq, brewer_name),
            ylab = ylab,
            xlab = xlab,
            ...)
    axis(side = 4, at = info, labels = names(info)) # info on quantiles
    par(mar = oldpar, las = 0)
}

