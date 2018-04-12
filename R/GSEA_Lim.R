#' ES for 2 tail GSEA (Lim et al., 2009) - Regulons
#'
#' This computes one overall enrichment score for two complementary gene sets on a given signature
#' as described in Lim et al., Pac. Symp. Biocomput. 2009, 504-515.
#'
#' @param signature Named vector of gene-level statistics (preferably t-statistics)
#' @param regulon list object as returned by aracne2regulon (viper package)
#' @param weight Weight for each gene-level statistic (default = 1)
#' @return An overall enrichment score for two complementary gene sets on a signature
#' @export
gsea2ES_regulon <- function(signature, regulon, weight = 1) {
    tfm <- regulon[['tfmode']]
    pos_targ <- names(tfm)[tfm>=0]
    neg_targ <- names(tfm)[tfm<0]
    position_pos <- match(pos_targ, names(sort(signature, decreasing = TRUE))) # 1st decreasing for pos. targets
    position_neg <- match(neg_targ, names(sort(signature, decreasing = FALSE))) # 2nd increasing for neg. targets
    names(position_pos) <- pos_targ
    names(position_neg) <- neg_targ
    rlistr <- rank(signature)
    rlistr <- rlistr[!(names(signature) %in% names(tfm))] # ranks of non-hits
    rlistr <- sort(c(rlistr, position_pos, position_neg))
    rlist <- signature[match(names(rlistr), names(signature))]
    x <- which(names(rlist) %in% names(tfm))
    nr <- sum(abs(rlist[x])^weight)
    nh <- length(rlist)-length(x)
    es <- rep(-(1/nh),length(rlist))
    es[x] <- abs(rlist[x])^weight/nr
    rs <- cumsum(es)
    return(es = max(rs) + min(rs))
}

#' ES for 2 tail GSEA (Lim et al., 2009) - two character vectors
#'
#' This computes one overall enrichment score for two complementary gene sets on a given signature
#' as described in Lim et al., Pac. Symp. Biocomput. 2009, 504-515.
#'
#' @param signature Named vector of gene-level statistics (preferably t-statistics)
#' @param set1 character vector, i.e. the first gene set
#' @param set2 character vector, i.e. the second gene set
#' @param weight Weight for each gene-level statistic (default = 1)
#' @return An overall enrichment score for two complementary gene sets on a signature
#' @export
gsea2ES <- function(signature, set1, set2, weight = 1) {

    position1 <- match(set1, names(sort(signature, decreasing = TRUE))) # 1st decreasing for pos. targets
    position2 <- match(set2, names(sort(signature, decreasing = FALSE))) # 2nd increasing for neg. targets
    names(position1) <- set1
    names(position2) <- set2
    rlistr <- rank(signature)
    rlistr <- rlistr[!(names(signature) %in% c(set1, set2))] # ranks of non-hits
    rlistr <- sort(c(rlistr, position1, position2))
    rlist <- signature[match(names(rlistr), names(signature))]
    x <- which(names(rlist) %in% c(set1, set2))
    nr <- sum(abs(rlist[x])^weight) # normalizing factor for hits
    nh <- length(rlist)-length(x) # normalizing factor for non-hits
    es <- rep(-(1/nh),length(rlist))
    es[x] <- abs(rlist[x])^weight/nr
    rs <- cumsum(es)
    return(es = max(rs) + min(rs))
}

