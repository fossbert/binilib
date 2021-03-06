# Functions for normalization of gene expression data, specifically RNA-Seq count data

#' tpm
#'
#' This function computes tpm for human exonic read counts
#'
#' @param counts Numeric matrix of \strong{human} raw read counts per gene
#' @param log Logical, whether the data should be log2 transformed
#' @param offset Number to offset zero values when applying log transformation
#' @param type Character to indicate the gene identifier, defaults to gene symbol
#' @param species Character to indicate the genes' species origin, defaults to human
#' @return Numeric matrix of TPM data
#' @export
tpm <- function(counts, log = TRUE, offset = 0.25,
                type = c('Symbol', 'Entrez'),
                species = c('Human', 'Mouse'),
                filterZero = TRUE) {

    if (filterZero) d1 <- counts[rowSums(counts) > 0, ] else d1 <- counts
    type <- match.arg(type, choices = c('Symbol', 'Entrez'))
    species <- match.arg(species, choices = c('Human', 'Mouse'))

    # retrieve information on gene lengths, gene symbols and entrez ids
    data("geneLengths")
    glen <- geneLengths[[species]][[type]]

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
#' @param tpm logical, adjust for gene length before rank transformation
#' @param pmOne logical, whether to rescale values so they are distributed between -1 and +1
#' @param gausstrans logical, whether to rescale the values to quantiles of the gaussian distribution
#' @param ... further arguments passed to the tpm function
#' @return basic rank normalized matrix.
#' @export

basic_rank_norm <- function(counts,
                            tpm = FALSE,
                            pmOne = TRUE,
                            gausstrans = FALSE,
                            ...) {

    if (tpm) d1 <- tpm(counts, log = FALSE, ...) else d1 <- counts

    if(pmOne){
        # rank transform and distributes between -1 and 1
        # first for columns
        d1 <- apply(d1, 2, rank)/(nrow(d1) + 1)*2-1
        # then for genes
        d2 <- t(apply(d1, 1, rank)/(ncol(d1) + 1)*2-1)

        # possibly transform
        if(gausstrans) d2 <- qnorm(d2/2+.5)

    } else {
        # rank transform each sample/column
        d1 <- t(t(apply(d1, 2, rank, na.last="keep"))/(colSums(!is.na(counts))+1))
        # rank transform each gene/row
        d2 <- t(apply(d1, 1, rank, na.last="keep"))/(rowSums(!is.na(d1))+1)
    }
    rownames(d2) <- rownames(d1)
    return(d2)
}

#' Rescale a numeric vector to range between 0 and 1
#'
#' @param x a numeric vector
#' @param na_rm logical, whether to remove NA values in call to range (TRUE)
#' @return rescaled numeric vector (minimum is 0, maximum is 1)
#' @export

rescale01 <- function(x, na_rm = TRUE){

        x <- as.numeric(x)

        rng <- range(x, na.rm = na_rm)
        (x - rng[1])/(rng[2] - rng[1])
}


#' basic_signature
#'
#' This function prepares a gene expression matrix before metaVIPER. After rank
#' normalization of each sample, the row medians are subtracted and every column
#' is divided by the row median absolute deviation (MAD). It will filter out genes with a mad
#' of zero.
#'
#' @param eset Numeric matrix of gene expression (raw counts or pre-processed)
#' @param tpm logical, whether tpm normalization should be carried out (use if eset is in raw counts)
#' @param ... further arguments to tpm function
#' @return Numeric matrix with a basic signature for each sample
#' @export
#'
basic_signature <- function(eset, tpm = FALSE, ...){

    if(tpm) eset <- tpm(eset, ...)

    rmad <- apply(eset, 1, mad)
    keep <- rmad > 0

    rank <- apply(eset[keep, ], 2, rank)
    rmed <- apply(eset[keep, ], 1, median)

    return(((rank - rmed)/rmad[keep]))
}




#' Variance stabilization transformation for RNAseq data
#'
#' This function stabilizes the variance, transform the data and add shot noise to the data
#'
#' @param emat matrix containing the raw counts, with genes in rows and samples in columns
#' @param fitType Character string indicating the type of fit for the dispersion (see DESeq::estimateDispersions)
#' @param seed Integer indicating the fixed seed for random numbers, 0 for not setting the seed
#' @param noise Logical, indicating whether to add a shot of random noise to the data (splitting ties)
#' @return Expression matrix
#' @export

DEtransform <- function(emat,
                        fitType=c("parametric", "local"),
                        noise=TRUE,
                        seed=42) {

    # Safety nets
    if (any(is.na(emat)))
        stop("NA values are not allowed in the count matrix")
    if (any(emat < 0)) {
        stop("some values in matrix are negative")
    }

    # Check whether input are integers
    if (!is.integer(emat)) {
        if (any(round(emat) != emat)) {
            message("Rounding non-integer read counts to integers")
            emat <- round(emat)
        }
        message("converting counts to integer mode")
        mode(emat) <- "integer"
    }

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = emat,
                                          colData = data.frame(condition = rep(1, ncol(emat)),
                                                               row.names = colnames(emat)),
                                          design = ~ 1)

    if (seed>0) set.seed(seed)
    fitType <- match.arg(fitType)
    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- DESeq2::estimateDispersions(dds, fitType=fitType)
    vmat <- DESeq2::getVarianceStabilizedData(dds)
    tmp <- vmat

    if (noise) {
        tmp <- unlist(apply(vmat, 2, function(i) {
            x <- sort(unique(i))
            x <- cbind(x[1:(length(x)-1)], x[2:length(x)])
            x <- cbind(x[, 1], sqrt(binilib::rowVar(x)))
            return(list(x))
        }), recursive=FALSE)

        tmp <- cbind(unlist(lapply(tmp, function(i) i[, 1]), use.names=F), unlist(lapply(tmp, function(i) i[, 2]), use.names=F))
        tmp1 <- stats::smooth.spline(tmp[, 1], tmp[, 2], spar=.5)
        tmp[tmp[, 1]>tmp1$x[which.min(tmp1$y)], 2] <- 0
        tmp1 <- stats::smooth.spline(tmp[, 1], tmp[, 2], spar=.5)
        tmp <- vmat+rnorm(length(vmat))*stats::predict(tmp1, vmat)$y
    }
    if (!is.null(colnames(tmp))) colnames(tmp) <- colnames(vmat)
    return(tmp)
}
