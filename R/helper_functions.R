#' Extract rows or columns from a matrix with order of interest
#'
#' @param inputmat matrix with row and column names
#' @param reflist character vector of row or column names to be extracted
#' @param margin where to look for the reflist, by default rows
#' @return new matrix with rows or columns in order of the reference list
extract_features <- function(inputmat, reflist, margin = 1) {

    mat <- as.matrix(inputmat)

    if (margin == 1) {
        tmp <- which(reflist %in% rownames(mat))
        available_names <- reflist[tmp]
    } else {
        tmp <- which(reflist %in% colnames(mat))
        available_names <- reflist[tmp]
    }

    if (length(available_names) == 0) {
        stop("No items found. Check whether names in reflist
             are in agreement with feature names")
    }

    message(length(available_names)," out of ",length(reflist)," names could be found!")

    if (margin == 1) return(mat[available_names, ])
    else return(mat[ ,available_names])

}


#' Collapsing many-to-one relationships for gene expression matrices
#'
#' @param inputmat matrix whose rows can partly be summarized using a factor
#' @param fac factor used for summarizing rows
#' @param method method of summary, defaults to mean
#' @param verbose whether to inform user regarding progress
#' @return new matrix with rows or columns in order of the reference list

collapse_multi <- function(inputmat, fac, method = c("mean", "median", "sum"), verbose = TRUE) {

    method <- match.arg(method)
    # ensure appropriate data types
    mat <- as.matrix(mat)
    fac <- as.factor(fac)

    outmat <- matrix(nrow = length(levels(fac)), ncol = ncol(mat),
                     dimnames = list(levels(fac), colnames(mat)))

    # message methods of collapsing, number of samples and date
    if (verbose) {
        message("\nCollapsing mutiple probes/IDs using their ", method, " for ", ncol(mat), " samples.")
        message("\nProcess started at ", date())
        pb <- txtProgressBar(max=ncol(mat), style = 3)
    }

    # message methods of collapsing, number of samples and date
    for (i in 1:ncol(mat)) {
        switch(method,
               mean={outmat[,i] <- tapply(mat[,i], fac, mean)},
               median={outmat[,i] <- tapply(mat[,i], fac, median)},
               sum={outmat[,i] <- tapply(mat[,i], fac, sum)})
        if (is(pb, "txtProgressBar")) setTxtProgressBar(pb, i)
    }

    if (verbose) message("\nProcess ended at ", date(), "\n")
    return(outmat)
}



