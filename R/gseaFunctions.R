#' This function performs Gene Set Enrichment Analysis as described by Subramanian et al.
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
gsea_one <- function(signature,
                     gS,
                     weight = 1,
                     sorting = c('decreasing', 'increasing'),
                     onlyES = FALSE){

                # sorting type
                sort_type <- match.arg(arg = sorting, choices = c('decreasing', 'increasing'))
                switch(
                    decreasing = {signature <- sort(signature, decreasing = TRUE)},
                    increasing = {signature <- sort(signature)}
                    )

                idx <- which(names(signature) %in% gS) # positions
                Nr <- sum(abs(signature[idx])^weight) # normalization factor
                Nh <- length(signature) - length(idx) # non-hits
                tmp <- rep(-(1/Nh), length(signature))
                tmp[idx] <- abs(signature[idx])^weight/Nr
                rs <- cumsum(tmp) # running sum
                maxabs <- which.max(abs(rs))
                es <- rs[maxabs]

                if(onlyES) es

                else {

                # leading edge
                 pos1 <- which.max(abs(rs))
                    if (es > 0) {
                        leg <- names(signature)[1:pos1]
                        leg <- leg[leg %in% gS]
                        legidx <- which(names(signature) %in% leg)
                    }
                    else {
                        leg <- names(signature)[pos1:length(signature)]
                        leg <- leg[leg %in% gS]
                        legidx <- which(names(signature) %in% leg)
                    }

                # return object
                gsea1 <- list(ES = es,
                              RS = rs,
                              signature = signature,
                              es_idx = pos1,
                              gs_idx = idx,
                              ledge = leg,
                              ledge_index = legidx)
                class(gsea1) <- "gsea1"
                return(gsea1)
                }

    }

#' Calculating a null distribution of enrichment scores
#'
#' This function generates a GSEA null distribution of enrichment scores based on gene shuffling, i.e.
#' the positions of genes in the signature are permuted rather than phenotype labels.
#'
#' @param signature named numeric vector of gene-level statistics
#' @param gS character vector of gene names matching names from the signature
#' @param w numeric indicating the weights used for GSEA (defaults to 1)
#' @param perm integer indicating the number of permutations to carry out (defaults to 1000)
#' @param seed integer for random number generation during the permutation process
#' @param plot_hist logical on whether to plot a histogram of the permutation-based enrichment scores
#' @param verbose logical on whether to message progress report
#' @return A gsea_null object
#' @export
gsea_null <- function(signature,
                      gS,
                      w = 1,
                      perm = 1000,
                      seed = 42,
                      plot_hist = FALSE,
                      verbose = TRUE) {

        set.seed(seed)

        # set up
        null_list <- lapply(1:perm, function(i, signature){
                tmp <- signature
                names(tmp) <- sample(names(signature))
                tmp
            }, signature = signature)

        # message date and what's going on
        if (verbose) {
            message("\nCalculating null distribution of enrichment scores by gene shuffling. Started at ", date())
            pb <- txtProgressBar(max = length(null_list), style = 3)
        } else pb <- NULL

        null_es <- vapply(seq_along(null_list), function(i, gS, w){

            if (!is.null(pb)) setTxtProgressBar(pb, i)

            gsea_one(signature = null_list[[i]], gS = gS, weight = w, onlyES = TRUE)

        }, FUN.VALUE = numeric(1), gS = gS, w = w)

        if(plot_hist) graphics::hist(null_es, breaks = perm/10, main = '')

        gsea_null <- list(null_es = null_es,
                          pos_null_es = mean(null_es[null_es >= 0]),
                          neg_null_es = mean(null_es[null_es < 0]))

        class(gsea_null) <- "gsea_null"
        return(gsea_null)
}


#' This function performs GSEA for the negative and positive targets, respectively, of a regulatory gene.
#'
#' @param signature Named vector of gene-level statistics (preferably t-statistics)
#' @param regulon list object as returned by aracne2regulon (viper package)
#' \describe{
#' \item{tfmode}{named numeric vector of the mode of regulation (MOR) values}
#' \item{likelihood}{numeric vector of interaction confidence values (based on mutual information)}
#' }
#' @param w weight for each gene-level statistic (default = 1)
#' @return a list object containing for each gene set:
#' \describe{
#' \item{ES}{enrichment score}
#' \item{RS}{running sum}
#' \item{signature}{input signature}
#' \item{es_idx}{index of the enrichment score}
#' \item{gs_idx}{indices of the gene set members in the signature}
#' \item{ledge}{gene identifiers of the genes in the leading edge}
#' \item{ledge_idx}{indices of the leading edge gene set members}
#' }
#' @export
gsea_regulon <- function(signature, regulon, weight = 1){

    # first isolate the two parts of the regulon
    tfmode <- regulon[[1]]
    pos_targ <- names(tfmode[tfmode >= 0])
    neg_targ <- names(tfmode[tfmode < 0])

    # little safety net
    if(!any(names(tfmode) %in% names(signature))) stop('Gene set names could not be found in signature!')

    # run gsea for each part of the regulon
    pos_gsea <- gsea_one(signature, pos_targ, weight = weight)
    neg_gsea <- gsea_one(signature, neg_targ, weight = weight)

    gsea.obj <- list(
        signature = pos_gsea$signature, # original signature
        ES_pos = pos_gsea$ES, ES_neg = neg_gsea$ES, # enrichment scores
        ES_pos_idx = pos_gsea$es_idx, ES_neg_idx = neg_gsea$es_idx, # index of enrichment scores
        RS_pos = pos_gsea$RS, RS_neg = neg_gsea$RS, # running sum for each gene set
        gs_idx_pos = pos_gsea$gs_idx, gs_idx_neg = neg_gsea$gs_idx, # indices for gene sets in signature
        ledge_pos = pos_gsea$ledge, ledge_neg = neg_gsea$ledge, # leading edges
        ledge_idx_pos = pos_gsea$ledge_index, ledge_idx_neg = neg_gsea$ledge_index # leading edge positions in signature
)

    class(gsea.obj) <- "gsea2"
    return(gsea.obj)
}


#' Plot results from a GSEA using a set of negative and positive targets of a regulatory gene
#'
#' This function generates a plot
#'
#' @param gsea2 gsea2 object as returned from \code{gsea_regulon}
#' @param color Vector of two components indicating the colors for the negative and positive parts of the regulon
#' @param ... Given for compatibility to the plot generic function
#' @return Nothing, a plot is generated in the default output device
#' @method plot gsea2
#' @export
plot.gsea2 <- function(gsea2, color=c("cornflowerblue","salmon"), ...) {

    # getting acquainted
    ylims <- (abs(round(max(c(gsea2$ES_pos, gsea2$ES_neg)), 1)) + 0.2)*c(-1, 1)

    # negative targets
    boxlim <-  c(max(ylims)-0.1, max(ylims))*sign(gsea2$ES_neg)
    plot(gsea2$RS_neg, col = color[1], las = 1, lwd = 2,
         ylim = ylims,
         type = 'l',
         xaxt = 'n',
         xlab = 'Gene signature index',
         ylab = 'ES',
         main = "",
         ...)
    abline(0, 0)
    for(i in 1:length(gsea2$gs_idx_neg)){
        lines(x = rep(seq_along(gsea2$signature)[gsea2$gs_idx_neg[i]], 2), y = boxlim, col = color[1], lwd = .5)
    }
    rect(xleft = 0, ybottom = boxlim[1], xright = length(gsea2$signature), ytop = boxlim[2])

    # positive targets
    boxlim <- c(max(ylims)-0.1, max(ylims))*sign(gsea2$ES_pos)
    lines(gsea2$RS_pos, col = color[2], las = 1, lwd = 2)
    for(i in 1:length(gsea2$gs_idx_pos)){
        lines(x = rep(seq_along(gsea2$signature)[gsea2$gs_idx_pos[i]], 2), y = boxlim, col = color[2], lwd = .5)
    }
    rect(xleft = 0, ybottom = boxlim[1], xright = length(gsea2$signature), ytop = boxlim[2])

}
