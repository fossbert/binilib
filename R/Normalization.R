# Functions for normalization of gene expression data, specifically RNA-Seq count data

#' tpm
#'
#' This function computes tpm for human exonic read counts
#'
#' @param counts Numeric matrix of raw read counts per gene
#' @param log Logical, whether the data should be log2 transformed
#' @param offset Number to offset zero values when applying log transformation
#' @return Numeric matrix of TPM data
#' @export
tpm <- function(counts, log = TRUE, offset = 0.25, type = c('Symbol', 'Entrez'), filterZero = TRUE) {

    if (filterZero) d1 <- counts[rowSums(counts) > 0, ] else d1 <- counts
    type <- match.arg(type, choices = c('Symbol', 'Entrez'))

    # retrieve information on gene lengths, gene symbols and entrez ids
    data("geneLengths")
    glen <- geneLengths[[type]]

    # match up
    common <- intersect(names(glen), rownames(d1))
    if (length(common) == 0) stop('No overlap of gene names could be established')
    d1 <- d1[common, ]
    glenKB <- glen[common]/1000

    # divide counts per samples by gene length in kb
    d2 <- d1/glenKB
    # get scaling factors
    scaling <- colSums(d2)/1e6
    # divide each column by its scaling factor
    d2 <- sweep(d2, 2, scaling, "/")

    if (!log) return(d2) else return(log2(d2 + offset))

}


#' Basic rank normalization of a gene expression matrix
#'
#' @param counts Matrix containing raw read counts (generally).
#' @param gausstrans whether or not the ranks should be transformed using a gaussian.
#' @return basic rank normalized matrix.
#' @export

basic_rank_norm <- function(counts, tpm = TRUE) {

    if (tpm) d1 <- tpm(counts, log = FALSE) else d1 <- counts

    # rank transform each sample/column
    d1 <- t(t(apply(d1, 2, rank, na.last="keep"))/(colSums(!is.na(counts))+1))
    rm(counts)
    gc()
    # rank transform each gene/row
    d2 <- t(apply(d1, 1, rank, na.last="keep"))/(rowSums(!is.na(d1))+1)
    rownames(d2) <- rownames(d1)
    return(d2)
    # rank transform and distributes between -1 and 1
    # first for columns
    #d1 <- apply(d1, 2, rank)/(nrow(d1) + 1)*2-1
    # then for genes
    #d2 <- t(apply(d1, 1, rank)/(ncol(d1) + 1)*2-1)

    # possibly transform
    # if (gausstrans) d2 <- qnorm(d2/2+.5)

}



