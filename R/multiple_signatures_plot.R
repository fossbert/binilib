#' Plot two tails of one regulon on multiple signatures
#'
#' This function generates a plot
#'
#' @param sigmat a matrix with genes in rows and signatures in columns. Typically a Z-score, t-statistic of log2 fold change matrix
#' @param tf a character vector of a regulatory protein, used to subset regulon object
#' @param regulon named regulon object similar to output of viper::aracne2regulon
#' @param color Vector of two components indicating the colors for the negative and positive parts of the regulon
#' @param nes numeric vector of normalized enrichment scores, must match number of columns in sigmat
#' @param padj numeric vector of adjusted p-values, must match number of columns in sigmat
#' @return Nothing, a plot is generated in the default output device
#' @export

plot_OneReg_MultSig <- function(sigmat,
                          tf,
                          regulon,
                          color=c("cornflowerblue", "salmon"),
                          nes = NULL,
                          padj = NULL,
                          ...) {

    omar <- par()$mar
    omgp <- par()$mgp

    tmp <- apply(sigmat, 2, gsea_regulon, regulon = regulon[[tf]])

    x <- seq(nrow(sigmat))
    y <- ncol(sigmat)

    layout(cbind(1,2,3), widths = c(1,1,6))

    # Plot 1: NES and FDR matrix
    #if(is.null)
    if(is.null(nes)) nes <- vapply(tmp, function(i) i$nes, FUN.VALUE = numeric(1))

    # change order for plot
    ordx <- order(nes)
    nes <- nes[ordx]

    mx <- max(nes, na.rm = TRUE)
    if(mx > 5) brks <- seq(-mx, mx, length.out = 101) else brks <- seq(-5, 5, length.out = 101)
    par(mar = c(2.1, 2.1, 4.1, .05))
    hcols <- c('royalblue', 'white', 'firebrick2')
    image(1, seq(y), t(nes), ylab="", xlab="",
          col = colorRampPalette(hcols)(length(brks)-1),
          breaks = brks,
          axes = FALSE, yaxs = "i")
    text(rep(1, y), seq(y), round(nes, 2), col = 'black', font = 2)
    box()
    grid(1, y, col="black", lty=1)
    mtext('NES', at = 1, line = 1)

    # Plot 2: FDR
    if(is.null(padj)) padj <- p.adjust(vapply(tmp, function(i) i$pval, FUN.VALUE = numeric(1)), 'fdr')
    padj <- padj[ordx]
    par(mar = c(2.1, 0.05, 4.1, 1.1))
    plot(NA, xlim = c(0, 1), ylim = c(0, y),
         type = 'n', axes = FALSE, ylab = "", xlab = "", yaxs = "i")
    text(rep(.5, y), seq(y)-.5, signif(padj, 2), col = 'black')
    grid(1, y, col="black", lty=1)
    mtext('FDR', at = .5, line = 1)

    # Plot 3: targets in signature
    xlimit <- c(0, length(x)*(1+.15*max(nchar(colnames(sigmat)))/8))
    tmp <- tmp[ordx]

    par(mar = c(2.1, 1.1, 4.1, 1.1))
    plot(0, type="n",
         ylim=c(0, y),
         xlim = xlimit,
         axes=FALSE,
         ylab="",
         xlab="",
         yaxs="i")
    # bar plots
    for (i in seq_along(tmp)) {

        gs2 <- tmp[[i]]
        blim_up <- c(i-.5, i)
        blim_down <- c(i-1, i-.5)

        xpos <- gs2$gs_idx_pos
        xneg <- gs2$gs_idx_neg

        for(j in seq_along(xpos)){
            lines(x = rep(x[xpos[j]], 2), y = blim_up, col = color[2], lwd = .8)
        }
        for(k in seq_along(xneg)){
            lines(x = rep(x[xneg[k]], 2), y = blim_down, col = color[1], lwd = .8)
        }
        abline(v = length(x))
        grid(0, y, col = 'black')

        #rect(xleft = 0, ybottom = min(c(blim_down, blim_up)),
         #    xright = length(x), ytop = max(c(blim_down, blim_up)))

    }
    text(rep(length(x)*1.02, y), seq(y)-.5, colnames(sigmat)[ordx], adj = 0)
    mtext(paste(tf, 'regulon'), at = floor(length(x)/2), line = 1)
    mtext('Signature', at = length(x), adj = 0, line = 1)

    par(mar = omar, mfrow = c(1,1), mfcol = c(1,1), las = 0)
}


#' Plot one-tailed GSEA on one signature for multiple gene sets
#'
#' This function generates a plot of the member genes of a given gene set on a given signature
#' for multiple gene sets.
#'
#' @param ges gene expression signature, a named vector of Z-scores, t-statistics, or log2 fold changes
#' @param geneSets a named list of character vectors representing the gene sets
#' @param fgseaRes Output of fgsea, will be used to get NES and adjusted p values
#' @param ledge_only Logical, whether to include all genes of gene set or only those in leading edge
#' @param signatureNames Character vector of length 2, specifying the experimental conditions
#' @param color Vector of two components indicating the colors for each part of the signature
#' @param ... adjustment of cex for pathway names
#' @return Nothing, a plot is generated in the default output device
#' @export

plot_fgseaRes <- function(ges,
                          geneSets,
                          fgseaRes,
                          ledge_only = FALSE,
                          signatureNames = NULL,
                           color=c("firebrick2", "royalblue"),
                           cex_pw = 1) {

    omar <- par()$mar
    omgp <- par()$mgp

    tmp <- lapply(geneSets, function(i){
        gsea1T(signature = ges, gS = i)
    })

    x <- seq(length(ges))
    y <- length(tmp)

    layout(cbind(1,2,3), widths = c(1,1,6))

    #    layout(rbind(c(1,1,2), c(3,4,5)),
 #          widths = c(1,1,6),
  #         heights = c(1, 5))

  #  if(plotSignature){
   #     layout(rbind(1,2),heights=c(1,3))
    #    par(mar = c(0, 3.1, 0.1, 1), mgp = c(2, .7, 0))
    #    plot(x, y,
    #         bty = 'n',
    #         type = "n",
    #         xaxt = 'n',
    #         yaxs = 'i',
    #         xlab = '',
    #         ylab = signatureType,
    #         las = 1,
    #         cex.axis = 0.8,
    #         cex.lab = 0.8)
    #
     #   switch(sigsort,
      #         '-1' = {
      #             polygon(c(min(x), x[y <= 0]), c(0, y[y <= 0]), border = NA, col = "gray80")
      #             polygon(c(max(x), x[y > 0]), c(0, y[y > 0]), border = NA, col = "gray80")
      #             if(!is.null(signatureNames)) {
      #                 text(min(x), max(y)/5, label = signatureNames[1], pos = 4)
      #                 text(max(x), min(y)/5, label = signatureNames[2], pos = 2)
      #             }
      #         },
      #         '1' = {
      #             polygon(c(min(x), x[y > 0]), c(0, y[y > 0]), border = NA, col = "gray80")
      #             polygon(c(max(x), x[y <= 0]), c(0, y[y <= 0]), border = NA, col = "gray80")
      #             if(!is.null(signatureNames)) {
      #                 text(min(x), min(y)/5, label = signatureNames[1], pos = 4)
      #                 text(max(x), max(y)/5, label = signatureNames[2], pos = 2)
      #             }
      #         })
      #  abline(h = 0, lty = 2, lwd = 2)
#    }

    # Get statistics from fgsea output

    if(!any(names(geneSets) %in% fgseaRes$pathway)){
        stop('Could not find any results for indicated gene sets!')
    }

    common <- intersect(names(tmp), fgseaRes$pathway)

    nes <- fgseaRes$NES
    names(nes) <- fgseaRes$pathway
    nes <- nes[common]
    padj <- fgseaRes$padj
    names(padj) <- fgseaRes$pathway
    padj <- padj[common]
    tmp <- tmp[common]

    # change order for plot
    ordx <- order(nes)
    nes <- nes[ordx]
    padj <- padj[ordx]
    tmp <- tmp[ordx]

    mx <- max(nes, na.rm = TRUE)
    if(mx > 5) brks <- seq(-mx, mx, length.out = 101) else brks <- seq(-5, 5, length.out = 101)
    par(mar = c(2.1, 2.1, 4.1, .05))
    hcols <- c('royalblue', 'white', 'firebrick2')
    image(1, seq(y), t(nes), ylab="", xlab="",
          col = colorRampPalette(hcols)(length(brks)-1),
          breaks = brks,
          axes = FALSE, yaxs = "i")
    text(rep(1, y), seq(y), round(nes, 2), col = 'black', font = 2)
    box()
    grid(1, y, col="black", lty=1)
    mtext('NES', at = 1, line = 1)

    # Plot 2: FDR
    par(mar = c(2.1, 0.05, 4.1, 1.1))
    plot(NA, xlim = c(0, 1), ylim = c(0, y),
         type = 'n', axes = FALSE, ylab = "", xlab = "", yaxs = "i")
    text(rep(.5, y), seq(y)-.5, signif(padj, 2), col = 'black')
    grid(1, y, col="black", lty=1)
    mtext('FDR', at = .5, line = 1)

    # Plot 3: targets in signature
    xlimit <- c(0, length(x)*(1+.15*max(nchar(names(tmp)))/8))

    par(mar = c(2.1, 0.1, 4.1, 1.1))
    plot(0, type="n",
         ylim=c(0, y),
         xlim = xlimit,
         axes=FALSE,
         ylab="",
         xlab="",
         yaxs="i")
    # bar plots

    # set up colors
    cols <- t(apply(col2rgb(color), 1, function(i){
        seq(i[1], i[2], length.out = length(x))
    }))

    for (i in seq_along(tmp)) {

        # tune alphas depending on NES sign

        if(sign(nes[i]) == 1){
            alphas <- seq(255, 50, length.out = length(x))
        } else{
            alphas <- seq(50, 255, length.out = length(x))
        }

        gs1 <- tmp[[i]]
        if(ledge_only) {
            gene_pos <- gs1$ledge_index
        } else gene_pos <- gs1$gs_idx
        blim <- c(i-1, i)

        for(j in seq_along(gene_pos)){
            idx <- gene_pos[j]
            tmpcol <- cols[ ,idx]
            barcol <- rgb(red = tmpcol[1],
                          green = tmpcol[2],
                          blue = tmpcol[3],
                          maxColorValue = 255,
                          alpha = alphas[idx])
            lines(x = rep(x[idx], 2), y = blim,
                  col = barcol, lwd = .8)
        }
        abline(v = length(x))
        grid(0, y, col = 'black')

        #rect(xleft = 0, ybottom = min(c(blim_down, blim_up)),
        #    xright = length(x), ytop = max(c(blim_down, blim_up)))

    }
    text(rep(length(x)*1.02, y), seq(y)-.5, names(tmp), adj = 0, cex = cex_pw)

    if(!is.null(signatureNames)){
        if(length(signatureNames) != 2) stop('Need 2 (!) signature names!')
        mtext(signatureNames[1], at = 0, line = 1, adj = 0)
        mtext(signatureNames[2], at = length(x), line = 1, adj = 1)
    }
    mtext('Pathway', at = (length(x) + xlimit[2])/2, line = 1)

    par(mar = omar, mfrow = c(1,1), mfcol = c(1,1), las = 0)
}


