# Vectorized convenience functions for matrix operations

#' rowVar
#'
#' This function replaces the apply(A, 1, var) call
#'
#' @param A a numeric matrix
#' @param removeNA should NA values be removed from rowSums/rowMeans ?
#' @return vector of variance per row
#' @export
rowVar <- function(A, removeNA = TRUE) {
    rowSums((A - rowMeans(A, na.rm = removeNA))^2, na.rm = removeNA)/(dim(A)[2] - 1)
}

#' colVar
#'
#' This function replaces the apply(A, 2, var) call
#'
#' @param A a numeric matrix
#' @param removeNA should NA values be removed from colSums/colMeans ?
#' @return vector of variance per column
#' @export
colVar <- function(A, removeNA = TRUE) {
    colSums((A - colMeans(A, na.rm = removeNA))^2, na.rm = removeNA)/(dim(A)[1] - 1)
}
