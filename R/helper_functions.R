#' Extract rows or columns from a matrix in order of interest
#'
#' @param expmat matrix with row and column names
#' @param reflist character vector of row or column names to be extracted
#' @param margin where to look for the reflist, by default rows
#' @return new matrix with rows or columns in order of the reference list
#' @export
extract_features <- function(expmat, reflist, margin = 1) {

    mat <- as.matrix(expmat)

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
#' @param expmat matrix whose rows can partly be summarized using a factor
#' @param fac factor used for summarizing rows
#' @param method method of summary, defaults to mean
#' @param verbose whether to inform user regarding progress
#' @return new matrix with rows or columns in order of the reference list
#' @export

collapse_multi <- function(expmat, fac, method = c("mean", "median", "sum"), verbose = TRUE) {

    method <- match.arg(method)
    # ensure appropriate data types
    mat <- as.matrix(expmat)
    fac <- as.factor(fac)

    outmat <- matrix(nrow = length(levels(fac)), ncol = ncol(mat),
                     dimnames = list(levels(fac), colnames(mat)))

    # message methods of collapsing, number of samples and date
    if (verbose) {
        message("\nCollapsing mutiple probes/IDs using their ", method, " for ", ncol(mat), " samples.")
        message("\nProcess started at ", date())
        pb <- txtProgressBar(max = ncol(mat), style = 3)
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

#' Finding the most variable features of a matrix
#'
#' @param expmat numerical matrix, typically a gene expression matrix
#' @param n number of most variable features to be returned, default 2000
#' @param method method to determine 'variable', defaults to robust measure (mad)
#' @param center whether to median center the resulting matrix, defaults to FALSE
#' @return new matrix with most variable genes
#' @export

find_var_genes <- function(expmat, n = 2000, method = c("mad","sd"), center = FALSE) {

    mat <- as.matrix(expmat)
    method <- match.arg(method)
    switch(method,
           mad={spread <- apply(mat, 1, mad)},
           sd={spread <- apply(mat, 1, sd)}
           )

    vmat <- mat[order(spread, decreasing = TRUE)[1:n], ]
    if (isTRUE(center)) {
        vmat <- sweep(vmat, 1, apply(vmat, 1, median, na.rm = TRUE))
    }
    return(vmat)
}


#' Tidyverse way to read in gene expression matrices
#'
#' @param df path to a tab-separated dataframe whose first column includes gene identifiers
#' @return Data matrix ready to play
#' @export
read_matrix <- function(df) {
        tmp <- readr::read_tsv(df)
        emat <- data.matrix(tmp[-1])
        rownames(emat) <- tmp[[1]]
        return(emat)
}

#' Convenience function for writing gene expression to file
#'
#' @param eset Numeric gene expression matrix with gene identifiers as row names
#' @param path character string indicating path and name of file
#' @export
write_matrix <- function(eset, path) {

    tmp <- tibble::rownames_to_column(data.frame(eset), var = 'Gene')

    readr::write_tsv(tmp, path = path)
}

#' lapply function to apply functions in a pair-wise manner between all distinct elements of a list
#'
#' @param x list
#' @param f the function to be applied between each pair of elements of the list. Function needs to take two arguments: i and j.
#' @return list object with the results
#' @export
#' @examples
#' lapply_pair(x = l, f = function(i, j) diag(cor(i, j)))
#' This call will return the column-wise Pearson correlation between each pair

lapply_pair <- function(x, f, ...) {

    if (is.null(names(x))) nom <- paste('Element', 1:length(x), sep = '_') else nom <- names(x)
    out <- vector("list", length(x))
    names(out) <- nom

    for (i in seq_along(x)) {

        for(j in seq_along(x)) {
            if (i == j) next
            id <- names(x)[j]
            out[[i]][[id]] <- f(x[[i]], x[[j]], ...)
        }

    }
    out
}


#' randomize_eset
#'
#' randomize columns or rows, respectively of a matrix
#'
#' @param eset Numeric matrix with row and column names
#' @param seed integer, defaults to 1
#' @param margin where to randomize, defaults to columns
#' @param ... Further arguments passed to apply
#' @return new matrix with randomized columns or rows, respectively
#' @export

randomize_eset <- function(eset, seed = 1, margin = 2, ...) {
         set.seed(1)

        if (margin == 2) {

            if (!is.null(rownames(eset))) nom <- rownames(eset)
            tmp <- apply(eset, margin, function(i, ...) sample(i, size = length(i), ...), ...)
            if (!is.null(rownames(eset))) rownames(tmp) <- nom

        } else {

            if (!is.null(colnames(eset))) nom <- colnames(eset)
            tmp <- t(apply(eset, margin, function(i, ...) sample(i, size = length(i), ...), ...))
            if (!is.null(colnames(eset))) colnames(tmp) <- nom
        }
       tmp
    }


#' Replace missing values (typically in a vector or matrix)
#'
#' @param x Vector (or matrix) with NA values somewhere
#' @param replacement Defaults to zero for numeric data, but can be anything really
#' @return new matrix with randomized columns or rows, respectively
#' @export

replace_missings <- function(x, replacement = 0) {

    if(!identical(class(x), class(replacement))) {
        warning('Classes of input and replacement do not match!')
    }
    is_miss <- is.na(x)
    x[is_miss] <- replacement
    message(sum(is_miss), ' missing value(s) replaced by the value ', replacement)
    x

}

#' PCA and with covariates
#'
#' Quickly explore whether categorical covariates are associated with major sources
#' of variance in the expression data.
#'
#' @param expmat matrix with row and column names
#' @param reflist character vector of row or column names to be extracted
#' @param margin where to look for the reflist, by default rows
#' @return new matrix with rows or columns in order of the reference list
#' @export
