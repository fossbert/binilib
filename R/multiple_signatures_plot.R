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
                          color=c("salmon", "cornflowerblue"),
                          nes = NULL,
                          padj = NULL,
                          ...) {

    omar <- par()$mar
    omgp <- par()$mgp

    tmp <- apply(sigmat, 2, gsea_regulon, regulon = regulon[[tf]])

    x <- seq(nrow(sigmat))
    y <- ncol(sigmat)

    # x-axis intervals
    xax_itvl <- .get_xintervals(x)

    layout(cbind(1,2,3), widths = c(1,1,6))

    # Plot 1: NES and FDR matrix

    if(is.null(nes)) nes <- vapply(tmp, function(i) i$nes, FUN.VALUE = numeric(1))

    # change order for plot
    ordx <- order(nes)
    nes <- nes[ordx]

    mx <- max(abs(nes), na.rm = TRUE)
    if(mx > 5) brks <- seq(-mx, mx, length.out = 101) else brks <- seq(-5, 5, length.out = 101)
    par(mar = c(2.1, 1.1, 2.1, .05))
    hcols <- c(color[2], 'white', color[1])
    image(1, seq(y), t(nes), ylab="", xlab="",
          col = colorRampPalette(hcols)(length(brks)-1),
          breaks = brks,
          axes = FALSE, yaxs = "i")
    text(rep(1, y), seq(y), round(nes, 2), col = 'black', ...)
    box()
    grid(1, y, col="black", lty=1)
    mtext('NES', at = 1, line = 0, ...)

    # Plot 2: FDR
    if(is.null(padj)) padj <- p.adjust(vapply(tmp, function(i) i$pval, FUN.VALUE = numeric(1)), 'fdr')
    padj <- padj[ordx]
    par(mar = c(2.1, 0.05, 2.1, 0.1))
    plot(NA, xlim = c(0, 1), ylim = c(0, y),
         type = 'n', axes = FALSE, ylab = "", xlab = "", yaxs = "i")
    text(rep(.5, y), seq(y)-.5, signif(padj, 2), col = 'black', ...)
    grid(1, y, col="black", lty=1)
    mtext('FDR', at = .5, line = 0, ...)

    # Plot 3: targets in signature
    xlimit <- c(0, length(x)*(1+.15*max(nchar(colnames(sigmat)))/8))
    tmp <- tmp[ordx]

    par(mar = c(2.1, 0.1, 2.1, 1.1))
    plot(0, type="n",
         ylim=c(0, y),
         xlim = xlimit,
         axes=FALSE,
         ylab="",
         xlab="",
         yaxs="i")
    # bar plots

    # set up barcolors
    cols <- col2rgb(color)

    for (i in seq_along(tmp)) {

        # tune alphas to blunt the middle
        middle <- floor(length(x)/2)
        alphas <- numeric(length(x))
        alphas[seq(middle)] <- seq(255, 100, length.out = middle)
        alphas[(middle+1):length(alphas)] <- seq(100, 255, length.out = length(alphas)-middle)

        gs2 <- tmp[[i]]
        blim_up <- c(i-.5, i-.2)
        blim_down <- c(i-.8, i-.5)

        xpos <- gs2$gs_idx_pos
        xneg <- gs2$gs_idx_neg

        tmpcol <- cols[,1]
        for(j in seq_along(xpos)){
            idxpos <- xpos[j]
            barcol <- rgb(red = tmpcol[1],
                          green = tmpcol[2],
                          blue = tmpcol[3],
                          maxColorValue = 255,
                          alpha = alphas[idxpos])
            lines(x = rep(x[idxpos], 2),
                  y = blim_up, col = barcol, lwd = .8)
        }

        tmpcol <- cols[,2]
        for(k in seq_along(xneg)){
            idxneg <- xneg[k]
            barcol <- rgb(red = tmpcol[1],
                          green = tmpcol[2],
                          blue = tmpcol[3],
                          maxColorValue = 255,
                          alpha = alphas[idxneg])
            lines(x = rep(x[idxneg], 2), y = blim_down, col = barcol, lwd = .8)
        }

        lines(x = c(0, length(x)), y = rep(blim_down[2], 2))
        rect(xleft = 0,
             ybottom = blim_down[1],
             xright = length(x),
             ytop = blim_up[2])

        # abline(v = length(x))
        # grid(0, y, col = 'black')

    }
    axis(side = 3, at = xax_itvl, tck = -.01, mgp=c(3,0,0),
         cex.axis = 1/2, ...)
    text(rep(length(x)*1.02, y), seq(y)-.5, colnames(sigmat)[ordx], adj = 0, ...)
    mtext(paste(tf, 'regulon'), side = 1, at = floor(length(x)/2), line = 0, ...)

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
#' @param strip_PWnames Logical, whether to remove commonly used Pathway prefixes, e.g. HALLMARK
#' @param ... adjustment of cex for pathway names
#' @return Nothing, a plot is generated in the default output device
#' @export

plot_fgseaRes <- function(ges,
                          geneSets,
                          fgseaRes,
                          ledge_only = FALSE,
                          signatureNames = NULL,
                           color=c("firebrick2", "royalblue"),
                          strip_PWnames = TRUE,
                          ...) {

    omar <- par()$mar
    omgp <- par()$mgp

    tmp <- lapply(geneSets, function(i){
        gsea1T(signature = ges, gS = i)
    })

    x <- seq(length(ges))

    # x-axis intervals
    xax_itvl <- .get_xintervals(x)

    layout(cbind(1,2,3), widths = c(1,1,6))

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
    y <- length(tmp)

    # change order for plot
    ordx <- order(nes)
    nes <- nes[ordx]
    padj <- padj[ordx]
    tmp <- tmp[ordx]

    mx <- max(abs(nes), na.rm = TRUE)
    if(mx > 3) brks <- seq(-mx, mx, length.out = 101) else brks <- seq(-3, 3, length.out = 101)
    par(mar = c(2.1, 1.1, 2.1, .05))
    hcols <- c(color[2], 'white', color[1])
    image(1, seq(y), t(nes), ylab="", xlab="",
          col = colorRampPalette(hcols)(length(brks)-1),
          breaks = brks,
          axes = FALSE,
          yaxs = "i")
    text(rep(1, y), seq(y), round(nes, 2), col = 'black', ...)
    box()
    grid(1, y, col="black", lty=1)
    mtext('NES', at = 1, line = 0, ...)

    # Plot 2: FDR
    par(mar = c(2.1, 0.05, 2.1, .05))
    plot(NA, xlim = c(0, 1), ylim = c(0, y),
         type = 'n', axes = FALSE, ylab = "", xlab = "", yaxs = "i")
    text(rep(.5, y), seq(y)-.5, signif(padj, 2), col = 'black', ...)
    grid(1, y, col="black", lty=1)
    mtext('FDR', at = .5, line = 0, ...)

    # Plot 3: targets in signature
    xlimit <- c(0, length(x)*(1+.2*max(nchar(names(tmp)))/8))

    par(mar = c(2.1, 0.05, 2.1, 1.1))
    plot(0,
         type="n",
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
        blim <- c(i-.8, i-.2)

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

        rect(xleft = 0,
             ybottom = blim[1],
             xright = length(x),
             ytop = blim[2])

    }
    if(strip_PWnames){
        pws <- gsub(pattern = '^[A-Z]+_', replacement = '', names(tmp))
    } else {
        pws <- names(tmp)
    }
    text(rep(length(x)*1.02, y), seq(y)-.5, pws, adj = 0, ...)
    axis(side = 3, at = xax_itvl, tck = -0.01, mgp=c(3,0,0),
         cex.axis = 1/2, ...)

    if(!is.null(signatureNames)){
      mtext(text = signatureNames[2], side = 1, line = 0, at = length(x), adj = 1, ...)
      mtext(text = signatureNames[1], side = 1, line = 0, at = 0, adj = 0, ...)
    }

    par(mar = omar, mfrow = c(1,1), mfcol = c(1,1), las = 0)
}




#' Plot ranks of a given gene set across several gene expression profiles from a numeric matrix.
#'
#' This function generates a plot of the member genes of a given gene set across multiple signatures.
#' Its first intention is to illustrate the so called 'percent rank' of a given gene from a given
#' gene set.
#'
#' @param sigmat gene expression matrix, usually of raw values such as read counts, TPM or FPKM
#' @param geneset a character vector representing one gene set
#' @param color single character string indicating a color for the bar lines
#' @param signatureNames Character vector of length 2, specifying the extremes of the signature
#' @param gs_label single character string to label the gene set under investigation
#' @param column_col vector of color ids (character, hex, etc) that must match the number of columns in the signature matrix.
#' @param scale_alpha logical, whether to scale alpha values for bars, defaults to TRUE
#' @param relative_sigs logical, whether signatures are considered relative, e.g. t-stats or logFC. Will carry out fgsea.
#' @param aREA_scores numeric vector that must match the number of columns from signature matrix, alternative to average rank metric
#' @param ... given for compatibility, particularly for the adjustment of cex for various text labels
#' @return Nothing, a plot is generated in the default output device
#' @examples
#' sigmat <- matrix(rnorm(1e5), ncol = 10, nrow = 1e4)
#' rownames(sigmat) <- paste0('Gene', seq(nrow(sigmat)))
#' colnames(sigmat) <- paste0('Sample', seq(ncol(sigmat)))
#' set.seed(42)
#' gs <- sample(rownames(sigmat), size = 100)
#' plot_OneGs_MultSig(sigmat = sigmat, geneset = gs)
#' @export

plot_OneGs_MultSig <- function(sigmat,
                               geneset,
                               color='black',
                               signatureNames = NULL,
                               gs_label = NULL,
                               column_col = NULL,
                               scale_alphas = TRUE,
                               relative_sigs = FALSE,
                               aREA_scores = NULL,
                               ...) {

  omar <- par()$mar
  omgp <- par()$mgp

  # x- and y-axis
  x <- seq(nrow(sigmat))
  y <- ncol(sigmat)

  tmp <- .get_index(ges_mat = sigmat, geneset = geneset)

    # order statistics
  if(relative_sigs){

    fres <- apply(sigmat, 2, function(i){
      fgsea::fgsea(pathways = list(gs = geneset), stats = i)
    }) %>% bind_rows()

    nes <- fres$NES
    fdrs <- p.adjust(fres$pval, 'fdr')
    ordx <- order(nes)

  } else {

    if(!is.null(aREA_scores)){
      if (length(aREA_scores) != ncol(sigmat)) stop('aREA scores need to match column number!')
      ordx <- order(aREA_scores)
      avg_score <- aREA_scores

    } else {
      score <- apply(sigmat, 2, function(i) rank(i, na.last = 'keep')/(length(!is.na(i))+1))
      common <- intersect(rownames(sigmat), geneset)
      avg_score <- colMeans(score[common, ], na.rm = TRUE)
      ordx <- order(avg_score)
    }

  }

  # change of order
  tmp <- tmp[ordx]

  # x-axis intervals
  xax_itvl <- .get_xintervals(x)

  # Plot 1: Avg pct rank matrix OR NES/FDRs

  if(relative_sigs){

    layout(cbind(1,2,3), widths = c(1,1,6))

    # change order for plot
    nes <- nes[ordx]
    fdrs <- fdrs[ordx]

    mx <- max(abs(nes), na.rm = TRUE)
    if(mx > 3) brks <- seq(-mx, mx, length.out = 101) else brks <- seq(-3, 3, length.out = 101)

     par(mar = c(2.1, 1.1, 2.1, .05))
    # defaults to red blue scheme
    hcols <- c('cornflowerblue', 'white', 'firebrick2')
    image(1, seq(y), t(nes), ylab="", xlab="",
          col = colorRampPalette(hcols)(length(brks)-1),
          breaks = brks,
          axes = FALSE,
          yaxs = "i")
    text(rep(1, y), seq(y), round(nes, 2), col = 'black', ...)
    box()
    grid(1, y, col="black", lty=1)
    mtext('NES', at = 1, line = 0, ...)

    # Plot 2: FDR
    par(mar = c(2.1, 0.05, 2.1, .05))
    plot(NA, xlim = c(0, 1), ylim = c(0, y),
         type = 'n', axes = FALSE, ylab = "", xlab = "", yaxs = "i")
    text(rep(.5, y), seq(y)-.5, signif(fdrs, 2), col = 'black', ...)
    grid(1, y, col="black", lty=1)
    mtext('FDR', at = .5, line = 0, ...)

  } else {

    layout(cbind(1,2), widths = c(1,6))

    # change order for plot
    avg_score <- avg_score[ordx]
    par(mar = c(2.1, 1.1, 2.1, 0.1))
    plot(0,
         type="n",
         ylim=c(0, y),
         axes=FALSE,
         ylab="",
         xlab="",
         yaxs="i")
    text(rep(1, y), seq(y)-.5, round(avg_score, 2), ...)
    box()
    grid(1, y, col="black", lty=1)
    if(!is.null(aREA_scores)){
      mtext(text = 'NES', side = 2, line = 0, ...)
    } else{
      mtext(text = 'Average percent rank', side = 2, line = 0, ...)
    }
  }

  # Plot 3: targets in signature
  xlimit <- c(0, length(x)*(1+.15*max(nchar(colnames(sigmat)))/8))

  par(mar = c(2.1, 0.1, 2.1, 1.1))
  plot(0, type="n",
       ylim=c(0, y),
       xlim = xlimit,
       axes=FALSE,
       ylab="",
       xlab="",
       yaxs="i")

  # set up alphas to be smallest in the middle
  if(scale_alphas){
    ap1 <- rank(x, na.last = 'keep')/(length(!is.na(x)) + 1)
    ap2 <- abs(ap1 - 0.5) * 2
    ap3 <- cut(ap2, breaks = quantile(ap2, probs = seq(0, 1, .05)), labels = FALSE, include.lowest=TRUE)
    ap <- ap3*1/length(unique(ap3))
  } else ap <- rep(1, length(x)) # no scaling

  rgb_val <- col2rgb(col = color)[,1]/255

  # bar plots
  for (i in seq_along(tmp)) {

    xpos <- tmp[[i]]
    blim <- c(i-.8, i-.2)

    for(j in seq_along(xpos)){
      idxpos <- xpos[j]
      barcol <- rgb(rgb_val[1],
                    rgb_val[2],
                    rgb_val[3],
                    alpha = ap[idxpos])

      lines(x = rep(x[idxpos], 2),
            y = blim,
            col = barcol,
            lwd = .8)

    }

    rect(xleft = 0,
         ybottom = blim[1],
         xright = length(x),
         ytop = blim[2])

  }
  axis(side = 3, at = xax_itvl, tck = -0.01, mgp=c(3,0,0),
       cex.axis = 1/2, ...)

  if(is.null(column_col)) {
    text(rep(length(x)*1.02, y), seq(y)-.5, colnames(sigmat)[ordx], adj = 0, ...)
  } else {
    text(rep(length(x)*1.02, y), seq(y)-.5, colnames(sigmat)[ordx], col = column_col[ordx], adj = 0, ...)
  }

  if(!is.null(signatureNames)){
    mtext(text = signatureNames[2], side = 1, line = 0, at = length(x), adj = 1, ...)
    mtext(text = signatureNames[1], side = 1, line = 0, at = 0, adj = 0, ...)
  }

  if(!is.null(gs_label)) {
    mtext(text = paste0(gs_label,' (n=', length(geneset), ')'), at = floor(length(x)/2),
          side = 1, line = 0, ...)
  } else {
    mtext(text = paste0('Gene set (n=', length(geneset), ')'), at = floor(length(x)/2),
          side = 1, line = 0, ...)
  }
  par(mar = omar, mfrow = c(1,1), mfcol = c(1,1), las = 0)
}



.get_index <- function(ges_mat, geneset){
  lapply(seq(ncol(ges_mat)), function(i){
    si <- sort(ges_mat[,i], decreasing = TRUE)
    idx <- which(names(si) %in% geneset)
    return(idx)
    })
}



.get_xintervals <- function(ges){

  options <- c(1000, 2000, 4000)

  # test 1000, 2000 and 4000 intervals
  ratios <- floor(length(ges)/options)
  ratios[ratios<1] <- NA_real_
  itvl <- options[which.min(ratios)]

  xax_itvl <- c(0, seq(floor(length(ges)/itvl))*itvl)

  # if there is a substantial difference between the last interval and the
  # length of the signature, go all the way
  diff_end <- length(ges) - xax_itvl[length(xax_itvl)]
  if(diff_end > length(ges)*.1){
    xax_itvl[length(xax_itvl)] <- length(ges)
  }

  return(xax_itvl)
}

