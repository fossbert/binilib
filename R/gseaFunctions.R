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
                if(!any(gS %in% names(signature))) stop('Gene set names could not be found in signature!')

                # GSEA Subramanian style
                idx <- which(names(signature) %in% gS) # positions
                Nr <- sum(abs(signature[idx])^weight) # normalization factor
                Nh <- length(signature) - length(idx) # non-hits
                tmp <- rep(-(1/Nh), length(signature))
                tmp[idx] <- abs(signature[idx])^weight/Nr
                rs <- cumsum(tmp) # running sum
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

        # loop through permutations and calculate enrichment scores
        null_es <- vapply(seq_along(null_list), function(i, gS, w){

            if (!is.null(pb)) setTxtProgressBar(pb, i)

            gsea1T(signature = null_list[[i]], gS = gS, weight = w, onlyES = TRUE)

        }, FUN.VALUE = numeric(1), gS = gS, w = w)

        if(plot_hist) graphics::hist(null_es, breaks = perm/10, main = '', xlab = 'Enrichment Score')

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
    if(!any(names(tfmode) %in% names(signature))) stop('Gene set names could not be found in signature!')

    # sorting type
    sort_type <- match.arg(arg = sorting, choices = c('decreasing', 'increasing'))
    # run gsea for each part of the regulon
    pos_gsea <- gsea1T(signature, pos_targ, sorting = sort_type, weight = weight)
    neg_gsea <- gsea1T(signature, neg_targ, sorting = sort_type, weight = weight)

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
        ledge_idx_pos = set1_gsea$ledge_index, ledge_idx_neg = set2_gsea$ledge_index # leading edge positions in signature
    )

    class(gsea.obj) <- "gsea2"
    return(gsea.obj)
}

#' Plot results from a one-tailed GSEA
#'
#' This function generates a plot
#'
#' @param gsea1 gsea1 object as returned from \code{gsea1T}
#' @param color character or integer indicating color choice (defaults to cornflowerblue)
#' @param main Character vector given as argument for convenience
#' @param plotSignature Logical whether to include a graph on the signature (e.g. t-statistics)
#' @param signatureNames Character vector indicating the conditions being compared. Defaults to NULL. Otherwise a charater vector of length two.
#' @param signatureType Character to indicate the type of gene level statistic (e.g. t-statistic).
#' Defaults to NULL which will yield a 'signature score' label.
#' @param ... Given for compatibility to the plot generic function
#' @return Nothing, a plot is generated in the default output device
#' @method plot gsea1
#' @export
plot.gsea1 <- function(gsea1,
                       color = 'cornflowerblue',
                       main = '',
                       plotSignature = FALSE,
                       signatureNames = NULL,
                       signatureType = NULL,
                       ...) {

    omar <- par()$mar
    omgp <- par()$mgp

    # brief set up
    x <- seq_along(gsea1$signature)
    y <- gsea1$signature
    if(is.null(signatureType)) signatureType <- 'signature score'

    if(plotSignature){
        layout(rbind(1,2),heights=c(1,3))
        par(mar = c(0, 3.1, 0.1, 1), mgp = c(2, .7, 0))
        plot(x, y,
             bty = 'n',
             type = "n",
             xaxt = 'n',
             yaxs = 'i',
             xlab = '',
             ylab = signatureType,
             las = 1,
             cex.axis = 0.8,
             cex.lab = 0.8)
        sigsort <- as.character(sign(y[1]))
        switch(sigsort,
               '-1' = {
                   polygon(c(min(x), x[y <= 0]), c(0, y[y <= 0]), border = NA, col = "gray80") ;
                   polygon(c(max(x), x[y > 0]), c(0, y[y > 0]), border = NA, col = "gray80")
                   if(!is.null(signatureNames)) {
                       text(min(x), max(y)/5, label = signatureNames[1], pos = 4)
                       text(max(x), min(y)/5, label = signatureNames[2], pos = 2)
                       }
                   },
               '1' = {
                   polygon(c(min(x), x[y > 0]), c(0, y[y > 0]), border = NA, col = "gray80")
                   polygon(c(max(x), x[y <= 0]), c(0, y[y <= 0]), border = NA, col = "gray80")
                   if(!is.null(signatureNames)) {
                   text(min(x), min(y)/5, label = signatureNames[1], pos = 4)
                   text(max(x), max(y)/5, label = signatureNames[2], pos = 2)
                        }
                   })
        abline(h = 0, lty = 2, lwd = 2)
    }

    # getting acquainted
    essign <- as.character(sign(gsea1$ES))
    switch(essign,
           '-1' = {ylims <- c(gsea1$ES - 0.2, max(gsea1$RS))},
           '1' = {ylims <- c(min(gsea1$RS), gsea1$ES + 0.2)})

    switch(essign,
           '-1' = {boxlim <- c(min(ylims) + 0.1, min(ylims))},
           '1' = {boxlim <- c(max(ylims) - 0.1, max(ylims))})

    xax_itvl <- seq(floor(length(x)/2000))*2000
    xax_lbls <- if(!plotSignature && !is.null(signatureNames)) {
        c(signatureNames[1], xax_itvl[1:(length(xax_itvl)-1)], signatureNames[2])
    } else c(0, xax_itvl)

    # time to plot
    par(mar = c(3.1, 3.1, 1, 1), mgp = c(2, .7, 0), las = 1)
    plot(x, gsea1$RS,
         col = color, las = 1, lwd = 2,
         ylim = ylims,
         type = 'l',
         xaxt = 'n',
         xlab = 'Gene signature index',
         ylab = 'ES',
         main = main,
         cex.axis = 0.8,
         ...)
    axis(side = 1, at = c(1, xax_itvl), labels = xax_lbls, cex.axis = 0.7)
    abline(h = 0, lty = 2, lwd = 2)
    for(i in seq_along(gsea1$gs_idx)){
        lines(x = rep(x[gsea1$gs_idx[i]], 2), y = boxlim, col = color, lwd = .5)
    }
    rect(xleft = 0, ybottom = boxlim[1], xright = length(x), ytop = boxlim[2])
    par(mar = omar, mgp = omgp, mfrow = c(1,1), mfcol = c(1,1), las = 0) # restore old settings
}



#' Plot 2-tailed enrichment results
#'
#' This function generates a plot of the GSEA running sums of two gene sets on a given signature.
#'
#' @param gsea2 gsea2 object as returned from \code{gsea_regulon} or \code{gsea2T}
#' @param color Vector of two components indicating the colors for the negative and positive parts of the regulon
#' @param legend Logical on whether to draw a legend explaining the two parts of the regulon
#' @param main Character vector given as argument for convenience
#' @param plotSignature Logical whether to include a graph on the signature (e.g. t-statistics)
#' @param signatureNames Character vector indicating the conditions being compared. Defaults to NULL. Otherwise a charater vector of length two.
#' @param signatureType Character to indicate the type of gene level statistic (e.g. t-statistic).
#' Defaults to NULL which will yield a 'signature score' label.
#' @param ... Given for compatibility to the plot generic function
#' @return Nothing, a plot is generated in the default output device
#' @method plot gsea2
#' @export
plot.gsea2 <- function(gsea2,
                       color = c("cornflowerblue","salmon"),
                       legend = TRUE,
                       main = '',
                       plotSignature = FALSE,
                       signatureNames = NULL,
                       signatureType = NULL,
                       ...) {

    omar <- par()$mar
    omgp <- par()$mgp

    x <- seq_along(gsea2$signature)
    y <- gsea2$signature
    if(is.null(signatureType)) signatureType <- 'signature score'

    if (plotSignature){
        layout(rbind(1,2), heights = c(1,3))
        par(mar = c(0, 3.1, 0.25, 1), mgp = c(2, 0.7, 0))
        plot(x, y,
             bty = 'n',
             type = "n",
             xaxt = 'n',
             xlab = '',
             ylab = signatureType,
             yaxs = 'i',
             las = 1,
             cex.axis = 0.8,
             cex.lab = 0.8,
             tck = -0.05
             )
        sigsort <- as.character(sign(y[1]))
        switch(sigsort,
               '-1' = {
                   polygon(c(min(x), x[y <= 0]), c(0, y[y <= 0]), border = NA, col = "gray80") ;
                   polygon(c(max(x), x[y > 0]), c(0, y[y > 0]), border = NA, col = "gray80")
                   if(!is.null(signatureNames)) {
                       text(min(x), max(y)/5, label = signatureNames[1], pos = 4)
                       text(max(x), min(y)/5, label = signatureNames[2], pos = 2)
                   }
               },
               '1' = {
                   polygon(c(min(x), x[y > 0]), c(0, y[y > 0]), border = NA, col = "gray80")
                   polygon(c(max(x), x[y <= 0]), c(0, y[y <= 0]), border = NA, col = "gray80")
                   if(!is.null(signatureNames)) {
                       text(min(x), min(y)/5, label = signatureNames[1], pos = 4)
                       text(max(x), max(y)/5, label = signatureNames[2], pos = 2)
                   }
               })
        abline(h = 0, lty = 2, lwd = 2)
    }

    # getting acquainted
    ylims <- (round(max(abs(c(gsea2$ES_pos, gsea2$ES_neg))), 1) + 0.2) * c(-1, 1)
    xax_itvl <- seq(floor(length(x)/2000))*2000
    xax_lbls <- if(!plotSignature && !is.null(signatureNames)) {
        c(signatureNames[1], xax_itvl[1:(length(xax_itvl)-1)], signatureNames[2])
    } else c(0, xax_itvl)

    # negative targets
    boxlimNeg <-  c(max(ylims)-0.1, max(ylims)) * sign(gsea2$ES_neg)
    par(mar = c(3.1, 3.1, 1, 1), mgp = c(2, 0.7, 0), las = 1)
    plot(x, gsea2$RS_neg,
         col = color[1], las = 1, lwd = 2,
         ylim = ylims,
         type = 'l',
         xaxt = 'n',
         xlab = 'Gene signature index',
         ylab = 'ES',
         main = main,
         cex.axis = 0.8,
         tck = -.025,
         ...)
    axis(side = 1, at = c(1, xax_itvl), labels = xax_lbls, tck = -0.025, cex.axis = 0.7)
    abline(h = 0, lty = 2)
    for(i in seq_along(gsea2$gs_idx_neg)){
        lines(x = rep(x[gsea2$gs_idx_neg[i]], 2), y = boxlimNeg, col = color[1], lwd = .5)
    }
    rect(xleft = 0, ybottom = boxlimNeg[1], xright = length(x), ytop = boxlimNeg[2])

    # positive targets
    boxlimPos <- c(max(ylims)-0.1, max(ylims))*sign(gsea2$ES_pos)
    lines(x, gsea2$RS_pos, col = color[2], las = 1, lwd = 2)
    for(i in seq_along(gsea2$gs_idx_pos)){
        lines(x = rep(x[gsea2$gs_idx_pos[i]], 2), y = boxlimPos, col = color[2], lwd = .5)
    }
    rect(xleft = 0, ybottom = boxlimPos[1], xright = length(x), ytop = boxlimPos[2])

    if(legend){
        # find distances from the ES
        essign <- as.character(sign(gsea2$ES_pos))

        switch(essign,
               '1' = {where <- as.character(which.min(c(gsea2$ES_pos_idx, length(x) - gsea2$ES_neg_idx)))},
               '-1' = {where <- as.character(which.min(c(gsea2$ES_neg_idx, length(x) - gsea2$ES_pos_idx)))})

        switch(where,
               '1' = {
                legend(x = length(x), y = ylims[2]-.15, title = 'Targets',
                       legend = c('pos', 'neg'), col = rev(color),
                       lty = 1, horiz = TRUE, x.intersp = .7, y.intersp = .7, xjust = 1,
                          lwd = 2, seg.len = .7, cex = 0.8)
               },
               '2' = {
                   legend(x = 1, y = ylims[1]+.15, title = 'Targets', legend = c('pos', 'neg'),
                          col = rev(color),
                          lty = 1, horiz = TRUE, x.intersp = .7, y.intersp = .7, xjust = 0,
                          yjust = 0, lwd = 2, seg.len = .7, cex = 0.8)
               })
    }
    par(mar = omar, mgp = omgp, mfrow = c(1,1), mfcol = c(1,1), las = 0) # restore old settings
}


