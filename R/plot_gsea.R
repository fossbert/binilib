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
#' @param ledge logical whether to plot the leading edge line
#' @param fixY numeric vector specyfying y-axis limits, defaults to NULL
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
                       ledge = TRUE,
                       fixY = NULL,
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
                   polygon(c(min(x), x[y <= 0]), c(0, y[y <= 0]), border = NA, col = "gray80")
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

    # if ylims should be fixed
    if(!is.null(fixY)){
        ylims <- fixY
        if(min(ylims) == 0){
            boxlim <- c(max(ylims) - 0.1, max(ylims))
            essign <- "1"
        } else {
            boxlim <- c(min(ylims) + 0.1, min(ylims))
            essign <- "-1"
        }
    } else{
        # if ylims roam freely
        essign <- as.character(sign(gsea1$ES))
        switch(essign,
               '-1' = {
                   ylims <- c(gsea1$ES - 0.2, max(gsea1$RS))
                   boxlim <- c(min(ylims) + 0.1, min(ylims))
               },
               '1' = {
                   ylims <- c(min(gsea1$RS), gsea1$ES + 0.2)
                   boxlim <- c(max(ylims) - 0.1, max(ylims))
               })
    }

    xax_itvl <- seq(floor(length(x)/2000))*2000
    xax_lbls <- if(!plotSignature && !is.null(signatureNames)) {
        c(signatureNames[1], xax_itvl[1:(length(xax_itvl)-1)], signatureNames[2])
    } else c(0, xax_itvl)

    # time to plot
    par(mar = c(3.1, 3.1, 1, 1), mgp = c(2, .7, 0), las = 1)
    plot(x, gsea1$RS,
         col = color, las = 1, lwd = 3,
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
        lines(x = rep(x[gsea1$gs_idx[i]], 2), y = boxlim, col = color, lwd = .8)
    }
    rect(xleft = 0, ybottom = boxlim[1], xright = length(x), ytop = boxlim[2])

    if(ledge) lines(x = rep(gsea1$es_idx, 2), y = c(0, gsea1$ES), col = color, lty = 2) # leading edge line

    if(!is.null(gsea1$NES)) {
        switch(essign,
               '-1' = {
                   text(x = round(length(x)/2),
                        y = boxlim[1], pos = 3,
                        labels = paste0('NES = ', round(gsea1$NES, 2), ', p = ', signif(gsea1$pval, 3)))
               },
               '1' = {
                   text(x = round(length(x)/2),
                        y = boxlim[1], pos = 1,
                        labels = paste0('NES = ', round(gsea1$NES, 2), ', p = ', signif(gsea1$pval, 3)))
               })

    }
    par(mar = omar, mgp = omgp, mfrow = c(1,1), mfcol = c(1,1), las = 0) # restore old settings
}



#' Plot 2-tailed enrichment results
#'
#' This function generates a plot of the GSEA running sums of two gene sets on a given signature.
#'
#' @param gsea2 gsea2 object as returned from \code{gsea2T}
#' @param color Vector of two components indicating the colors for the negative and positive parts of the regulon
#' @param legend Logical on whether to draw a legend explaining the two parts of the regulon
#' @param legNames Character, defaults to Targets, pos and neg
#' @param main Character vector given as argument for convenience
#' @param plotSignature Logical whether to include a graph on the signature (e.g. t-statistics)
#' @param signatureNames Character vector indicating the conditions being compared. Defaults to NULL. Otherwise a charater vector of length two.
#' @param signatureType Character to indicate the type of gene level statistic (e.g. t-statistic).
#' @param ledge logical whether to plot the leading edge line
#' @param fixY numeric vector specyfying y-axis limits, defaults to NULL
#' Defaults to NULL which will yield a 'signature score' label.
#' @param ... Given for compatibility to the plot generic function
#' @return Nothing, a plot is generated in the default output device
#' @method plot gsea2
#' @export
plot.gsea2 <- function(gsea2,
                       color = c("cornflowerblue","salmon"),
                       legend = TRUE,
                       legNames = c('Targets', 'pos', 'neg'),
                       main = '',
                       plotSignature = FALSE,
                       signatureNames = NULL,
                       signatureType = NULL,
                       ledge = TRUE,
                       fixY = NULL,
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
    if(!is.null(fixY)){
        ylims <- fixY
    } else{
        ylims <- (round(max(abs(c(gsea2$ES_pos, gsea2$ES_neg))), 1) + 0.2) * c(-1, 1)
    }
    xax_itvl <- seq(floor(length(x)/2000))*2000
    xax_lbls <- if(!plotSignature && !is.null(signatureNames)) {
        c(signatureNames[1], xax_itvl[1:(length(xax_itvl)-1)], signatureNames[2])
    } else c(0, xax_itvl)

    # define boxes
    if(sign(gsea2$ES_pos) == sign(gsea2$ES_neg)) {
        # if ES have same sign, fix order for boxes, positive above, negative down
        boxlimPos <- c(max(ylims)-0.1, max(ylims))
        boxlimNeg <-  c(min(ylims)+0.1, min(ylims))

    } else {
        boxlimPos <- c(max(ylims)-0.1, max(ylims))*sign(gsea2$ES_pos)
        boxlimNeg <-  c(max(ylims)-0.1, max(ylims)) * sign(gsea2$ES_neg)
    }

    # time to plot - starting with negative targets
    par(mar = c(3.1, 3.1, 1, 1), mgp = c(2, 0.7, 0), las = 1)
    plot(x, gsea2$RS_neg,
         col = color[1], las = 1, lwd = 3,
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
    abline(h = 0, lty = 2, lwd = 2) # center line
    # bar tags for negative targets
    for(i in seq_along(gsea2$gs_idx_neg)){
        lines(x = rep(x[gsea2$gs_idx_neg[i]], 2), y = boxlimNeg, col = color[1], lwd = .8)
    }
    rect(xleft = 0, ybottom = boxlimNeg[1], xright = length(x), ytop = boxlimNeg[2])

    # running sum plot for positive targets
    lines(x, gsea2$RS_pos, col = color[2], las = 1, lwd = 3)

    # bar tags for positive targets
    for(i in seq_along(gsea2$gs_idx_pos)){
        lines(x = rep(x[gsea2$gs_idx_pos[i]], 2), y = boxlimPos, col = color[2], lwd = .8)
    }
    rect(xleft = 0, ybottom = boxlimPos[1], xright = length(x), ytop = boxlimPos[2])

    # leading edge lines
    if(ledge) {
        lines(x = rep(gsea2$ES_pos_idx, 2),
              y = c(0, gsea2$ES_pos), col = color[2],
              lty = 2)
        lines(x = rep(gsea2$ES_neg_idx, 2),
              y = c(0, gsea2$ES_neg), col = color[1],
              lty = 2)
    }

    # legend
    if(legend){
        # find distances from the ES
        essign <- as.character(sign(gsea2$ES_pos))

        switch(essign,
               '1' = {where <- as.character(which.min(c(gsea2$ES_pos_idx, length(x) - gsea2$ES_neg_idx)))},
               '-1' = {where <- as.character(which.min(c(gsea2$ES_neg_idx, length(x) - gsea2$ES_pos_idx)))})

        switch(where,
               '1' = {
                   legend(x = length(x), y = ylims[2]-.15, title = legNames[1],
                          legend = legNames[2:3], col = rev(color),
                          lty = 1, horiz = TRUE, x.intersp = .7, y.intersp = .7, xjust = 1,
                          lwd = 2, seg.len = .7, cex = 0.8)
               },
               '2' = {
                   legend(x = 1, y = ylims[1]+.15, title = legNames[1], legend = legNames[2:3],
                          col = rev(color),
                          lty = 1, horiz = TRUE, x.intersp = .7, y.intersp = .7, xjust = 0,
                          yjust = 0, lwd = 2, seg.len = .7, cex = 0.8)
               })
    }
    # NES and p-value if available
    if(!is.null(gsea2$NES_pos)) {
        # again in case NES have same sign ... fix order
        if(sign(gsea2$NES_pos) == sign(gsea2$NES_neg)){
            text(x = round(length(x)/2), y = boxlimPos[1], pos = 1,
                 labels = paste0('NES = ', round(gsea2$NES_pos, 2), ', p = ', signif(gsea2$pval_pos, 3)))
            text(x = round(length(x)/2), y = boxlimNeg[1], pos = 3,
                 labels = paste0('NES = ', round(gsea2$NES_neg, 2), ', p = ', signif(gsea2$pval_neg, 3)))

            # otherwise put NES according to ES direction
        } else {
            if(gsea2$NES_pos >= 0){
                text(x = round(length(x)/2), y = boxlimPos[1], pos = 1,
                     labels = paste0('NES = ', round(gsea2$NES_pos, 2), ', p = ', signif(gsea2$pval_pos, 3)))
                text(x = round(length(x)/2), y = boxlimNeg[1], pos = 3,
                     labels = paste0('NES = ', round(gsea2$NES_neg, 2), ', p = ', signif(gsea2$pval_neg, 3)))
            } else {
                text(x = round(length(x)/2), y = boxlimPos[1], pos = 3,
                     labels = paste0('NES = ', round(gsea2$NES_pos, 2), ', p = ', signif(gsea2$pval_pos, 3)))
                text(x = round(length(x)/2), y = boxlimNeg[1], pos = 1,
                     labels = paste0('NES = ', round(gsea2$NES_neg, 2), ', p = ', signif(gsea2$pval_neg, 3)))
            }
        }
    }
    par(mar = omar, mgp = omgp, mfrow = c(1,1), mfcol = c(1,1), las = 0) # restore old settings
}







#' Plot 2-tailed enrichment results
#'
#' This function generates a plot of the GSEA running sums of two gene sets on a given signature.
#'
#' @param gsea2regulon gsea2regulon object as returned from \code{gsea_regulon}
#' @param color Vector of two components indicating the colors for the negative and positive parts of the regulon
#' @param legend Logical on whether to draw a legend explaining the two parts of the regulon
#' @param legNames Character, defaults to Targets, pos and neg
#' @param main Character vector given as argument for convenience
#' @param plotSignature Logical whether to include a graph on the signature (e.g. t-statistics)
#' @param signatureNames Character vector indicating the conditions being compared. Defaults to NULL. Otherwise a charater vector of length two.
#' @param signatureType Character to indicate the type of gene level statistic (e.g. t-statistic).
#' @param ledge logical whether to plot the leading edge line
#' @param fixY numeric vector specyfying y-axis limits, defaults to NULL
#' @param ... Given for compatibility to the plot generic function
#' @return Nothing, a plot is generated in the default output device
#' @method plot gsea2regulon
#' @export
plot.gsea2regulon <- function(gsea2regulon,
                       color = c("cornflowerblue","salmon"),
                       legend = TRUE,
                       legNames = c('Targets', 'pos', 'neg'),
                       main = '',
                       plotSignature = FALSE,
                       signatureNames = NULL,
                       signatureType = NULL,
                       ledge = TRUE,
                       fixY = NULL,
                       ...) {

    omar <- par()$mar
    omgp <- par()$mgp

    x <- seq_along(gsea2regulon$signature)
    y <- gsea2regulon$signature
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
    if(!is.null(fixY)){
        ylims <- fixY
    } else{
        ylims <- (round(max(abs(c(gsea2regulon$ES_pos, gsea2regulon$ES_neg))), 1) + 0.2) * c(-1, 1)
    }
    xax_itvl <- seq(floor(length(x)/2000))*2000
    xax_lbls <- if(!plotSignature && !is.null(signatureNames)) {
        c(signatureNames[1], xax_itvl[1:(length(xax_itvl)-1)], signatureNames[2])
    } else c(0, xax_itvl)

    # define boxes
    if(sign(gsea2regulon$ES_pos) == sign(gsea2regulon$ES_neg)) {
        # if ES have same sign, fix order for boxes, positive above, negative down
        boxlimPos <- c(max(ylims)-0.1, max(ylims))
        boxlimNeg <-  c(min(ylims)+0.1, min(ylims))

    } else {
        boxlimPos <- c(max(ylims)-0.1, max(ylims))*sign(gsea2regulon$ES_pos)
        boxlimNeg <-  c(max(ylims)-0.1, max(ylims)) * sign(gsea2regulon$ES_neg)
    }

    # time to plot - starting with negative targets
    par(mar = c(3.1, 3.1, 1, 1), mgp = c(2, 0.7, 0), las = 1)
    plot(x, gsea2regulon$RS_neg,
         col = color[1], las = 1, lwd = 3,
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
    abline(h = 0, lty = 2, lwd = 2) # center line
    # bar tags for negative targets
    for(i in seq_along(gsea2regulon$gs_idx_neg)){
        lines(x = rep(x[gsea2regulon$gs_idx_neg[i]], 2), y = boxlimNeg, col = color[1], lwd = .8)
    }
    rect(xleft = 0, ybottom = boxlimNeg[1], xright = length(x), ytop = boxlimNeg[2])

    # running sum plot for positive targets
    lines(x, gsea2regulon$RS_pos, col = color[2], las = 1, lwd = 3)

    # bar tags for positive targets
    for(i in seq_along(gsea2regulon$gs_idx_pos)){
        lines(x = rep(x[gsea2regulon$gs_idx_pos[i]], 2), y = boxlimPos, col = color[2], lwd = .8)
    }
    rect(xleft = 0, ybottom = boxlimPos[1], xright = length(x), ytop = boxlimPos[2])

    # leading edge lines
    if(ledge) {
        lines(x = rep(gsea2regulon$ES_pos_idx, 2),
              y = c(0, gsea2regulon$ES_pos), col = color[2],
              lty = 2)
        lines(x = rep(gsea2regulon$ES_neg_idx, 2),
              y = c(0, gsea2regulon$ES_neg), col = color[1],
              lty = 2)
    }

    # legend
    if(legend){
        # find distances from the ES
        essign <- as.character(sign(gsea2regulon$ES_pos))

        switch(essign,
               '1' = {where <- as.character(which.min(c(gsea2regulon$ES_pos_idx, length(x) - gsea2regulon$ES_neg_idx)))},
               '-1' = {where <- as.character(which.min(c(gsea2regulon$ES_neg_idx, length(x) - gsea2regulon$ES_pos_idx)))})

        switch(where,
               '1' = {
                   legend(x = length(x), y = ylims[2]-.15, title = legNames[1],
                          legend = legNames[2:3], col = rev(color),
                          lty = 1, horiz = TRUE, x.intersp = .7, y.intersp = .7, xjust = 1,
                          lwd = 2, seg.len = .7, cex = 0.8)
               },
               '2' = {
                   legend(x = 1, y = ylims[1]+.15, title = legNames[1], legend = legNames[2:3],
                          col = rev(color),
                          lty = 1, horiz = TRUE, x.intersp = .7, y.intersp = .7, xjust = 0,
                          yjust = 0, lwd = 2, seg.len = .7, cex = 0.8)
               })
    }

      # NES and p-value
    if(sign(gsea2regulon$nes) == sign(gsea2regulon$ES_pos)){
        text(x = round(length(x)/2),
             y = boxlimPos[1], pos = ifelse(gsea2regulon$ES_pos > 0, 1, 3),
             labels = paste0('NES = ',
                             round(gsea2regulon$nes, 2),
                             ', p = ',
                             signif(gsea2regulon$pval, 3)))
    } else {
        text(x = round(length(x)/2),
             y = boxlimNeg[1], pos = ifelse(gsea2regulon$ES_neg > 0, 1, 3),
             labels = paste0('NES = ',
                             round(gsea2regulon$nes, 2),
                             ', p = ',
                             signif(gsea2regulon$pval, 3)))
        }

    par(mar = omar, mgp = omgp, mfrow = c(1,1), mfcol = c(1,1), las = 0) # restore old settings
}


