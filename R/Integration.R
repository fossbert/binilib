#' Integration with trapezoid method
#'
#' This function integrate over a numerical range using the trapezoid method
#'
#' @param x Numeric vector of x values
#' @param y Numeric vector of y values
#' @return Number
#' @export

integrateTZ <- function(x, y) {
    pos <- order(x)
    x <- x[pos]
    y <- y[pos]
    idx = 2:length(x)
    return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

#' Numerical integration of functions
#'
#' Integrates numerically a function over a range using the trapezoid method
#'
#' @param f Function of 1 variable (first argument)
#' @param xmin Number indicating the min x value
#' @param xmax Number indicating the max x value
#' @param steps Integer indicating the number of steps to evaluate
#' @param ... Additional arguments for \code{f}
#' @return Number
#' @export

integrateFunction <- function(f, xmin, xmax, steps=100, ...) {
    x <- seq(xmin, xmax, length=steps)
    y <- f(x, ...)
    integrateTZ(x, y)
}

#' Integration based on CDF (AOC)
#'
#' This function integrates a distribution of scores based on the area over the CDF curve
#'
#' @param x Numeric vector, matrix or list of vectors or matrixes
#' @param xlim Numeric vector of 2 elements indicating the range where to perform the integration
#' @details This function computes the area over the curve for the vector o columns of the matrix provided as input
#' @export

cdfInteg <- function(x, xlim=NULL) {
    if (is.null(xlim)) xlim <- range(unlist(x, use.names=FALSE), na.rm=TRUE)
    if (is.list(x)) return(sapply(x, cdfInteg, xlim=xlim))
    if (is.matrix(x)) return(apply(x, 2, cdfInteg, xlim=xlim))
    1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim)
}
