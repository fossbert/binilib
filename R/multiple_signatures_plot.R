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

