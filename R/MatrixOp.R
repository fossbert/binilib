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
#' @param removeNA logical should NA values be removed from colSums/colMeans ?
#' @return vector of variance per column
#' @export
colVar <- function(A, removeNA = TRUE) {
    colSums((A - colMeans(A, na.rm = removeNA))^2, na.rm = removeNA)/(dim(A)[1] - 1)
}


#' Calculate row-wise Welch’s t-statistic
#'
#' @param expmat numerical matrix, typically a gene expression matrix
#' @param classes integer vector indicating class assignment (0 or 1)
#' @return numeric vector with Welch's t-statistics for each row
#' @export

rowTstat <- function(X, cl){
    X0 <- X[, cl == 0]
    X1 <- X[, cl == 1]
    m0 <- rowMeans(X0)
    m1 <- rowMeans(X1)
    n0 <- sum(cl == 0)
    n1 <- sum(cl)
    sq <- function(x) x * x
    s0 <- rowSums(sq(X0 - m0))
    s0 <- s0 / (n0 * (n0 - 1))
    s1 <- rowSums(sq(X1 - m1))
    s1 <- s1 / (n1 * (n1 - 1))
    (m0 - m1) / sqrt(s0 + s1)
}


