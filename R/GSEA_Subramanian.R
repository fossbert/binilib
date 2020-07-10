### These collection of functions are meant to implement classic GSEA functionality as described by
### Subramanian et al., PNAS, 2005

#' This function performs one-tailed Gene Set Enrichment Analysis as described by Subramanian et al.
#'
#' @param signature Named vector of gene-level statistics (preferably t-statistics)
#' @param gS Gene set of interest
#' @param w Weight for each gene-level statistic (default = 1)
#' @param sorting character indicating how to sort the signature (default is decreasing)
#' @param onlyES logical, whether only the enrichment score should be returned
#' @return A list of class "gsea1" containing:
#' \describe{
#' \item{ES}{enrichment score}
#' \item{RS}{running sum}
#' \item{signature}{input signature}
#' \item{es_idx}{index of the maximum enrichment score}
#' \item{gs_idx}{indices of the genes in the gene set in the signature}
#' \item{ledge}{gene identifiers for the genes in the leading edge}
#' \item{ledge_index}{indices of the leading edge genes in the signature}
#' }
#' @export
gsea1T <- function(signature,
                     gS,
                     weight = 1,
                     sorting = c('decreasing', 'increasing'),
                     onlyES = FALSE){

                # sorting type
                sort_type <- match.arg(arg = sorting, choices = c('decreasing', 'increasing'))

                switch(sort_type,
                    decreasing = {signature <- sort(signature, decreasing = TRUE)},
                    increasing = {signature <- sort(signature)}
                    )

                # little safety net
                if(!any(gS %in% names(signature))) stop('Gene set names could not be found in signature!', call. = FALSE)

                # GSEA Subramanian style
                idx <- which(names(signature) %in% gS) # positions
                Nr <- sum(abs(signature[idx])^weight) # normalization factor
                Nh <- length(signature) - length(idx) # non-hits
                tmp <- rep(-(1/Nh), length(signature))
                tmp[idx] <- abs(signature[idx])^weight/Nr
                rs <- cumsum(tmp) # running sum

                # plotting cosmetics:
                # adjust first and last position if first or last position is a hit
                if(rs[1] != -(1/Nh)) rs[1] <- -(1/Nh)
                if(rs[length(rs)] != -(1/Nh)) rs[length(rs)] <- -(1/Nh)

                maxabs <- which.max(abs(rs))
                es <- rs[maxabs]

                if(onlyES) es # for null model calculation

                else {

                # leading edge
                    if (es > 0) {
                        leg <- names(signature)[1:maxabs]
                        leg <- leg[leg %in% gS]
                        legidx <- which(names(signature) %in% leg)
                    }
                    else {
                        leg <- names(signature)[maxabs:length(signature)]
                        leg <- leg[leg %in% gS]
                        legidx <- which(names(signature) %in% leg)
                    }

                # return object
                gsea1 <- list(ES = es,
                              RS = rs,
                              signature = signature,
                              es_idx = maxabs,
                              gs_idx = idx,
                              ledge = leg,
                              ledge_index = legidx,
                              NES = NULL,
                              pval = NULL,
                              null_es = NULL,
                              mean_nulles = NULL)
                class(gsea1) <- "gsea1"
                return(gsea1)
                }
    }

#' This function performs GSEA for the negative and positive targets, respectively, of a regulatory gene.
#'
#' @param signature Named vector of gene-level statistics (preferably t-statistics)
#' @param regulon list object as returned by aracne2regulon (viper package)
#' \describe{
#' \item{tfmode}{named numeric vector of the mode of regulation (MOR) values}
#' \item{likelihood}{numeric vector of interaction confidence values (based on mutual information)}
#' }
#' @param sorting character indicating how to sort the signature (default is decreasing). Passed to gsea1T.
#' @param w weight for each gene-level statistic (default = 1)
#' @return a list object containing for each gene set:
#' \describe{
#' \item{signature}{input signature}
#' \item{ES_pos}{enrichment score positive targets}
#' \item{ES_neg}{enrichment score negative targets}
#' \item{ES_pos_idx}{enrichment score positive targets index}
#' \item{ES_neg_idx}{enrichment score negative targets index}
#' \item{RS_pos}{running sum for positive targets}
#' \item{RS_neg}{running sum for negative targets}
#' \item{gs_idx_pos}{indices of the positive targets in the signature}
#' \item{gs_idx_neg}{indices of the negative targets in the signature}
#' \item{ledge_pos}{gene identifiers of the leading edge for the positive targets}
#' \item{ledge_neg}{gene identifiers of the leading edge for the negative targets}
#' \item{ledge_idx_pos}{indices of the leading edge gene set members - positive targets}
#' \item{ledge_idx_neg}{indices of the leading edge gene set members - negative targets}
#' }
#' @export
gsea_regulon <- function(signature,
                         regulon,
                         sorting = c('decreasing', 'increasing'),
                         weight = 1){

    # first isolate the two parts of the regulon
    tfmode <- regulon[[1]]
    pos_targ <- names(tfmode[tfmode >= 0])
    neg_targ <- names(tfmode[tfmode < 0])

    # little safety net
    if(!any(names(tfmode) %in% names(signature))) stop('Gene set names could not be found in signature!', call. = FALSE)

    # sorting type
    sort_type <- match.arg(arg = sorting, choices = c('decreasing', 'increasing'))
    # run gsea for each part of the regulon
    pos_gsea <- gsea1T(signature, pos_targ, sorting = sort_type, weight = weight)
    neg_gsea <- gsea1T(signature, neg_targ, sorting = sort_type, weight = weight)

     #calc stats with aREA right away
    stats <- aREA_single(ges = signature, regulon = regulon)

    gsea.obj <- list(
        signature = pos_gsea$signature, # original signature
        ES_pos = pos_gsea$ES, ES_neg = neg_gsea$ES, # enrichment scores
        ES_pos_idx = pos_gsea$es_idx, ES_neg_idx = neg_gsea$es_idx, # index of enrichment scores
        RS_pos = pos_gsea$RS, RS_neg = neg_gsea$RS, # running sum for each gene set
        gs_idx_pos = pos_gsea$gs_idx, gs_idx_neg = neg_gsea$gs_idx, # indices for gene sets in signature
        ledge_pos = pos_gsea$ledge, ledge_neg = neg_gsea$ledge, # leading edges
        ledge_idx_pos = pos_gsea$ledge_index, ledge_idx_neg = neg_gsea$ledge_index, # leading edge positions in signature
        nes = stats$nes, pval = stats$pval
)

    class(gsea.obj) <- "gsea2regulon"
    return(gsea.obj)
}


#' This function performs GSEA of two gene sets for a given signature
#'
#' Similar to \code{gsea_regulon} but without any information on Mode of Regulation
#'
#' @param signature Named vector of gene-level statistics (preferably t-statistics)
#' @param set1 character vector, i.e. the first gene set
#' @param set2 character vector, i.e. the second gene set
#' @param sorting character indicating how to sort the signature (default is decreasing). Passed to gsea1T.
#' @param w weight for each gene-level statistic (default = 1)
#' @return a list object containing for each gene set:
#' \describe{
#' \item{signature}{input signature}
#' \item{ES_pos}{enrichment score set1}
#' \item{ES_neg}{enrichment score set2}
#' \item{ES_pos_idx}{enrichment score set1 index}
#' \item{ES_neg_idx}{enrichment score set2 index}
#' \item{RS_pos}{running sum for set1}
#' \item{RS_neg}{running sum for set2}
#' \item{gs_idx_pos}{indices of the set1 genes in the signature}
#' \item{gs_idx_neg}{indices of the set2 genes in the signature}
#' \item{ledge_pos}{gene identifiers of the leading edge for set1}
#' \item{ledge_neg}{gene identifiers of the leading edge for set2}
#' \item{ledge_idx_pos}{indices of the leading edge gene set members - set1}
#' \item{ledge_idx_neg}{indices of the leading edge gene set members - set2}
#' }
#' @export
gsea2T <- function(signature,
                   set1,
                   set2,
                   sorting = c('decreasing', 'increasing'),
                   weight = 1){

    # little safety net
    if(!any(c(set1, set2) %in% names(signature))) stop('Gene set names could not be found in signature!')

    # sorting type
    sort_type <- match.arg(arg = sorting, choices = c('decreasing', 'increasing'))

    # run gsea for each of the gene sets
    set1_gsea <- gsea1T(signature, set1, sorting = sort_type, weight = weight)
    set2_gsea <- gsea1T(signature, set2, sorting = sort_type, weight = weight)

    # return time ! For compatibility, list names are similar to gsea_regulon, although "positive"
    # and "negative" don't make too much sense
    gsea.obj <- list(
        signature = set1_gsea$signature, # original signature
        ES_pos = set1_gsea$ES, ES_neg = set2_gsea$ES, # enrichment scores
        ES_pos_idx = set1_gsea$es_idx, ES_neg_idx = set2_gsea$es_idx, # index of enrichment scores
        RS_pos = set1_gsea$RS, RS_neg = set2_gsea$RS, # running sum for each gene set
        gs_idx_pos = set1_gsea$gs_idx, gs_idx_neg = set2_gsea$gs_idx, # indices for gene sets in signature
        ledge_pos = set1_gsea$ledge, ledge_neg = set2_gsea$ledge, # leading edges
        ledge_idx_pos = set1_gsea$ledge_index, ledge_idx_neg = set2_gsea$ledge_index, # leading edge positions in signature
        null_es_pos = NULL, null_es_neg = NULL, NES_pos = NULL, NES_neg = NULL,
        pval_pos = NULL, pval_neg = NULL
    )

    class(gsea.obj) <- "gsea2"
    return(gsea.obj)
}
