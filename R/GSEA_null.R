#' @include classes.R
NULL

#' GSEA null model
#'
#' This function generates a GSEA null distribution of enrichment scores based on gene shuffling, i.e.
#' the positions of genes in the signature are permuted rather than phenotype labels. It will add a
#' normalized enrichment score (NES) to the gsea object as described by Subramanian et al. and statistical
#' inference is carried out using one-tailed tests on the appropriate (positive/negative)
#' side of the null distribution. For two-tailed GSEA, the null model is calculated accordingly and
#' separate NES and p-values are provided for each gene set.
#'
#' @param gsea_obj a gsea object as returned by gsea1T or gsea2T/gsea_regulon
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
setMethod("gsea_null", c(gsea_obj="gsea1"), function(gsea_obj,
                                         w = 1,
                                         perm = 1000,
                                         seed = 42,
                                         analytical = FALSE,
                                         verbose = TRUE) {

    signature <- gsea_obj$signature
    gs_size <- length(gsea_obj$gs_idx)

    set.seed(seed)

    # set up permutations
    gs0 <- lapply(seq(perm), function(i){
        sample(names(signature), size = gs_size)
    })

    # message date and what's going on
    if (verbose) {
        message("\nCalculating null distribution of enrichment scores by gene shuffling.\nStarted at ", date())
        pb <- txtProgressBar(max = perm, style = 3)
    } else pb <- NULL

    # loop through permutations and calculate enrichment scores
    null_es <- vapply(seq(perm), function(i, signature, w){

        if (!is.null(pb)) setTxtProgressBar(pb, i)

        gsea1T(signature = signature, gS = gs0[[i]], weight = w, onlyES = TRUE)

    }, FUN.VALUE = numeric(1), signature = signature, w = w)

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
            gsea_obj$pval <- (sum(null_es > gsea_obj$ES)+1)/(perm+1)
        } else {
            gsea_obj$pval <- (sum(null_es < gsea_obj$ES)+1)/(perm+1)
        }
    }
    return(gsea_obj)
})


#' @rdname gsea_null-methods
#' @aliases gsea_null, gsea2-method
setMethod("gsea_null", c(gsea_obj="gsea2"), function(gsea_obj,
                                         w = 1,
                                         perm = 1000,
                                         seed = 42,
                                         analytical = FALSE,
                                         verbose = TRUE) {

    signature <- gsea_obj$signature
    size_pos <- length(gsea_obj$gs_idx_pos)
    size_neg <- length(gsea_obj$gs_idx_neg)
    # how was the signature sorted, decreasing or increasing ?
    sigsort <- as.character(sign(signature[1]))
    set.seed(seed)

    # set up permutations
    gs0_pos <- lapply(seq(perm), function(i){
        sample(names(signature), size = size_pos)
    })

    gs0_neg <- lapply(seq(perm), function(i){
        sample(names(signature), size = size_neg)
    })

    # message date and what's going on
    if (verbose) {
        message("\nCalculating null distribution of enrichment scores by gene shuffling.\nStarted at ", date())
        pb <- txtProgressBar(max = perm, style = 3)
    } else pb <- NULL

    # loop through permutations and calculate enrichment scores
    null_es_pos <- vector(mode = 'numeric', length = perm)
    null_es_neg <- vector(mode = 'numeric', length = perm)

    for(i in seq(perm)){

        if (!is.null(pb)) setTxtProgressBar(pb, i)

        null_es_pos[i] <- gsea1T(signature = signature, gS = gs0_pos[[i]], weight = w, onlyES = TRUE)
        null_es_neg[i] <- gsea1T(signature = signature, gS = gs0_neg[[i]], weight = w, onlyES = TRUE)
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
            gsea_obj$pval_pos <- (sum(null_es_pos > gsea_obj$ES_pos)+1)/(length(null_es_pos)+1)
        } else gsea_obj$pval_pos <- (sum(null_es_pos < gsea_obj$ES_pos)+1)/(length(null_es_pos)+1)

        if(gsea_obj$ES_neg >= 0) {
            gsea_obj$pval_neg <- (sum(null_es_neg > gsea_obj$ES_neg)+1)/(length(null_es_neg)+1)
        } else gsea_obj$pval_neg <- (sum(null_es_neg < gsea_obj$ES_neg)+1)/(length(null_es_neg)+1)
    }
    return(gsea_obj)
})

