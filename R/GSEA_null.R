#' @include classes.R
NULL

#' Calculating a null distribution of enrichment scores for two-tail GSEA
#'
#' This function generates a GSEA null distribution of enrichment scores based on gene shuffling, i.e.
#' the positions of genes in the signature are permuted rather than phenotype labels.
#'
#' @param signature named numeric vector of gene-level statistics
#' @param set1 character vector of gene names of positive targets
#' @param set2 character vector of gene names of negative targets
#' @param w numeric indicating the weights used for GSEA (defaults to 1)
#' @param perm integer indicating the number of permutations to carry out (defaults to 1000)
#' @param seed integer for random number generation during the permutation process
#' @param plot_hist logical on whether to plot a histogram of the permutation-based enrichment scores
#' @param verbose logical on whether to message progress report
#' @return A gsea_null object
#' @export
gsea2T_null <- function(signature,
                        set1,
                        set2,
                        w = 1,
                        perm = 1000,
                        seed = 42,
                        plot_hist = FALSE,
                        verbose = TRUE) {

    set.seed(seed)

    # set up permutations
    null_list <- lapply(1:perm, function(i, signature){
        tmp <- signature
        names(tmp) <- sample(names(signature))
        tmp
    }, signature = signature)

    # message date and what's going on
    if (verbose) {
        message("\nCalculating null distribution of enrichment scores by gene shuffling.\nStarted at ", date())
        pb <- txtProgressBar(max = length(null_list), style = 3)
    } else pb <- NULL

    # little safety net
    if(!any(c(set1, set2) %in% names(signature))) stop('Gene set names could not be found in signature!', call. = FALSE)

    # loop through permutations and calculate enrichment scores
    null_es <- vapply(seq_along(null_list), function(i, set1, set2, w){

        if (!is.null(pb)) setTxtProgressBar(pb, i)

        gsea2ES(signature = null_list[[i]], set1, set2, weight = w)

    }, FUN.VALUE = numeric(1), set1 = set1, set2 = set2, w = w)

    if(plot_hist) graphics::hist(null_es, breaks = perm/10, main = '', xlab = 'Enrichment Score')

    gsea_null <- list(null_es = null_es,
                      pos_null_es = mean(null_es[null_es >= 0]),
                      neg_null_es = mean(null_es[null_es < 0]))

    class(gsea_null) <- "gsea_null2"
    return(gsea_null)
}


#' GSEA null model
#'
#' This function generates a GSEA null distribution of enrichment scores based on gene shuffling, i.e.
#' the positions of genes in the signature are permuted rather than phenotype labels. It will add a
#' normalized enrichment score (NES) to the gsea object as described by Subramanian et al. and statistical
#' inference is carried out using one-tailed tests on the appropriate (positive/negative)
#' side of the null distribution. For two-tailed GSEA, the null model is calculated accordingly and
#' separate NES and p-values are provided for each gene set.
#'
#' @param gsea_obj
#' @param w numeric indicating the weights used for GSEA (defaults to 1)
#' @param perm integer indicating the number of permutations to carry out (defaults to 1000)
#' @param seed integer for random number generation during the permutation process (defaults to 42)
#' @param analytical logical whether to derive p-values analytically, i.e. from two-tailed test using a normal distribution (defaults to FALSE, i.e. inference according to Subramanian et al.)
#' @param ... Additional parameters added to keep compatibility
#' @return gsea_obj with information on nullmodel added in extra slots
#' @export
#' @docType methods
#' @rdname gsea_null-methods
setGeneric("gsea_null", function(gsea_obj, ...) standardGeneric("gsea_null"))

#' @rdname gsea_null-methods
#' @aliases gsea_null,gsea1-method
setMethod("gsea_null", "gsea1", function(gsea_obj,
                                         w = 1,
                                         perm = 1000,
                                         seed = 42,
                                         analytical = FALSE,
                                         verbose = TRUE) {

    signature <- gsea_obj$signature
    gS <- names(signature)[gsea_obj$gs_idx]

    set.seed(seed)

    # set up permutations
    null_list <- lapply(1:perm, function(i, signature){
        tmp <- signature
        names(tmp) <- sample(names(signature))
        tmp
    }, signature = signature)

    # message date and what's going on
    if (verbose) {
        message("\nCalculating null distribution of enrichment scores by gene shuffling.\nStarted at ", date())
        pb <- txtProgressBar(max = length(null_list), style = 3)
    } else pb <- NULL

    # little safety net
    if(!any(gS %in% names(signature))) stop('Gene set names could not be found in signature!', call. = FALSE)

    # loop through permutations and calculate enrichment scores
    null_es <- vapply(seq_along(null_list), function(i, gS, w){

        if (!is.null(pb)) setTxtProgressBar(pb, i)

        gsea1T(signature = null_list[[i]], gS = gS, weight = w, onlyES = TRUE)

    }, FUN.VALUE = numeric(1), gS = gS, w = w)

    # assign null distribution of enrichment scores
    gsea_obj$null_es <- null_es
    # derive NES depending on sign of enrichment score
    if(gsea_obj$ES >= 0) {
        gsea_obj$mean_nulles <- mean(null_es[null_es >= 0])
    } else {
        gsea_obj$mean_nulles <- mean(null_es[null_es < 0])
    }
    gsea_obj$NES <- round(gsea_obj$ES/abs(gsea_obj$mean_nulles), 2)

    # derive p-values: 1.) analytically, two-tailed
    if(analytical){
        gsea_obj$pval <- signif(pnorm(abs(gsea_obj$NES), lower.tail = FALSE)*2, 3)
    } else { # 2.) permutation based, one-tailed
        if(gsea_obj$ES >= 0) {
            pval <- sum(null_es > gsea_obj$ES)/length(null_es)
        } else {
            pval <- sum(null_es < gsea_obj$ES)/length(null_es)
        }
        # adjust if p-value is 0
        if(pval == 0) gsea_obj$pval <- 1/length(null_es) else gsea_obj$pval <- pval

    }
    return(gsea_obj)
})


#' @rdname gsea_null-methods
#' @aliases gsea_null, gsea2-method
setMethod("gsea_null", "gsea2", function(gsea_obj,
                                         w = 1,
                                         perm = 1000,
                                         seed = 42,
                                         analytical = FALSE,
                                         verbose = TRUE) {

    signature <- gsea_obj$signature
    set_pos <- names(signature)[gsea_obj$gs_idx_pos]
    set_neg <- names(signature)[gsea_obj$gs_idx_neg]
    # how was the signature sorted, decreasing or increasing ?
    sigsort <- as.character(sign(signature[1]))
    set.seed(seed)

    # set up permutations
    null_list <- lapply(1:perm, function(i, signature){
        tmp <- signature
        names(tmp) <- sample(names(signature))
        tmp
    }, signature = signature)

    # message date and what's going on
    if (verbose) {
        message("\nCalculating null distribution of enrichment scores by gene shuffling.\nStarted at ", date())
        pb <- txtProgressBar(max = length(null_list), style = 3)
    } else pb <- NULL

    # little safety net
    if(!any(c(set_pos, set_neg) %in% names(signature))) stop('Gene set names could not be found in signature!', call. = FALSE)

    # loop through permutations and calculate enrichment scores
    null_es_pos <- vector(mode = 'numeric', length = perm)
    null_es_neg <- vector(mode = 'numeric', length = perm)

    for(i in seq_along(null_list)){

        if (!is.null(pb)) setTxtProgressBar(pb, i)

        null_es_pos[i] <- gsea1T(signature = null_list[[i]], gS = set_pos, weight = w, onlyES = TRUE)
        null_es_neg[i] <- gsea1T(signature = null_list[[i]], gS = set_neg, weight = w, onlyES = TRUE)
    }

    gsea_obj$null_es_pos <- null_es_pos
    gsea_obj$null_es_neg <- null_es_neg

    # identify the right tail for each enrichment score

    if(gsea_obj$ES_pos >= 0) {
       gsea_obj$NES_pos <- gsea_obj$ES_pos/mean(gsea_obj$null_es_pos[gsea_obj$null_es_pos >= 0])
    } else gsea_obj$NES_pos <- gsea_obj$ES_pos/abs(mean(gsea_obj$null_es_pos[gsea_obj$null_es_pos < 0]))

    if(gsea_obj$ES_neg >= 0) {
        gsea_obj$NES_neg <- gsea_obj$ES_neg/mean(gsea_obj$null_es_neg[gsea_obj$null_es_neg >= 0])
    } else gsea_obj$NES_neg <- gsea_obj$ES_neg/abs(mean(gsea_obj$null_es_neg[gsea_obj$null_es_neg < 0]))

    # derive p-value, either analytically or permutation based (default)

    if(analytical){
        gsea_obj$pval_pos <- pnorm(abs(gsea_obj$NES_pos), lower.tail = FALSE)*2
        gsea_obj$pval_neg <- pnorm(abs(gsea_obj$NES_neg), lower.tail = FALSE)*2
    } else {
        # permutation again depends on signs of enrichment scores
        if(gsea_obj$ES_pos >= 0) {
            pval_pos <- sum(null_es_pos > gsea_obj$ES_pos)/length(null_es_pos)
        } else pval_pos <- sum(null_es_pos < gsea_obj$ES_pos)/length(null_es_pos)

        if(gsea_obj$ES_neg >= 0) {
            pval_neg <- sum(null_es_neg > gsea_obj$ES_neg)/length(null_es_neg)
        } else pval_neg <- sum(null_es_neg < gsea_obj$ES_neg)/length(null_es_neg)

        # adjust if pval is 0
        if(pval_pos == 0) gsea_obj$pval_pos <- 1/length(null_es_pos) else gsea_obj$pval_pos <- pval_pos
        if(pval_neg == 0) gsea_obj$pval_neg <- 1/length(null_es_neg) else gsea_obj$pval_neg <- pval_neg
    }

    return(gsea_obj)
})

