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

#' Extract the top up and down regulated MR per sample
#'
#' @param vipres numeric matrix with results from VIPER
#' @param nn integer of number of top MR to consider
#' @param direction character, which tail should be returned, defaults to both up and down
#' @return list object with top dysregulated MR per sample
#' @export

topMR <- function(vipres, nn = 25, direction = c('both', 'up', 'down')) {

    dir <- match.arg(direction)

    switch(dir,
           both = {idx <- c(1:nn, (nrow(vipres)-(nn-1)):nrow(vipres))
                    tfm <- rep(c(1,-1), each = nn) },
           up = {idx <- 1:nn
                    tfm <- rep(1, nn)},
           down = {idx <- (nrow(vipres)-(nn-1)):nrow(vipres)
                    tfm <- rep(-1, nn)}
           )

        tmp <- apply(vipres, 2, function(i) {
            names(tfm) <- rownames(vipres[order(i, decreasing = TRUE), ])[idx]
            list(tfmode = tfm, likelihood = unname(abs(tfm)/nn))
        })

        class(tmp) <- 'regulon'
        return(tmp)
}


