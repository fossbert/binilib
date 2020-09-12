

# Functions for cluster analysis
#' Simple one-tail rank based enrichment analysis sREA
#'
#' This function performs simple 1-tail rank based enrichment analysis
#'
#' @param signatures Numeric matrix of signatures
#' @param groups List containing the groups as vectors of sample names
#' @return Matrix of Normalized Enrichment Zcores
#' @export
sREA <- function(signatures, groups) {
    if (is.null(nrow(signatures))) signatures <- matrix(signatures, length(signatures), 1, dimnames=list(names(signatures), "sample1"))
    sig <- qnorm(apply(signatures, 2, rank)/(nrow(signatures)+1))
    gr <- sapply(groups, function(x, samp) {
        samp %in% x
    }, samp=rownames(sig))
    gr <- t(gr)
    nn <- rowSums(gr)
    gr <- gr/nn
    es <- gr %*% sig
    return(es*sqrt(nn))
}

#' Cluster membership reliability estimated by enrichment analysis
#'
#'This function etimates the cluster membership reliability using aREA
#'
#' @param cluster Vector of cluster memberships or list of cluster memberships
#' @param similarity Similarity matrix
#' @param xlim Optional vector of 2 components indicating the limits for computing AUC
#' @param method Character string indicating the mthod to compute reliability, either by element, by cluster or global
#' @return Reliability score for each element
#' @export
clusterReliability <- function(cluster, similarity, xlim=NULL, method=c("element", "cluster", "global")) {
    method <- match.arg(method)
    if (!is.list(cluster)) cluster <- list(cluster=cluster)
    switch(method,
           "element" = {
               res <- lapply(cluster, function(cluster, similarity) {
                   reg <- split(names(cluster), cluster)
                   tmp <- sREA(similarity, reg)
                   tmp <- lapply(1:nrow(tmp), function(i, tmp, reg) {
                       tmp[i, ][colnames(tmp) %in% reg[[which(names(reg)==(rownames(tmp)[i]))]]]
                   }, tmp=tmp, reg=reg)
                   res <- unlist(tmp, use.names=FALSE)
                   names(res) <- unlist(lapply(tmp, names), use.names=FALSE)
                   return(res[match(names(cluster), names(res))])
               }, similarity=similarity)
               if (length(res)==1) return(res[[1]])
               return(res)
           },
           "cluster" = {
               rel <- clusterReliability(cluster, similarity, method="element")
               if (!is.list(rel)) rel <- list(cluster=rel)
               if (is.null(xlim)) xlim <- range(unlist(rel, use.names=FALSE))
               res <- lapply(1:length(cluster), function(i, cluster, rel, xlim) {
                   tapply(rel[[i]], cluster[[i]], function(x, xlim) {
                       1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim)
                   }, xlim=xlim)
               }, cluster=cluster, rel=rel, xlim=xlim)
               if (length(res)==1) return(res[[1]])
               return(res)
           },
           "global" = {
               rel <- clusterReliability(cluster, similarity, method="element")
               if (!is.list(rel)) rel <- list(cluster=rel)
               if (is.null(xlim)) xlim <- range(unlist(rel, use.names=FALSE))
               res <- lapply(rel, function(x, xlim) {
                   1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim)
               }, xlim=xlim)
               if (length(res)==1) return(res[[1]])
               return(res)
           }
    )
}

#' Plot cluster results
#'
#' This function generates a plot for the cluster results in the default output device
#'
#' @param dset Numeric matrix to include as heatmap
#' @param hdist Distance object for the columns
#' @param vdist Distance object for the rows
#' @param hcluster Vector with cluster classes for the columns
#' @param vcluster Vector with cluster classes for the rows
#' @param hclass Vector or matrix of horizontal classes
#' @param vclass Vector or matrix of vertical classes
#' @param hscore Optional vector of scores for the columns
#' @param vscore Optional vector of scores for the rows
#' @param hmargins Vector of 4 elements for the heatmap margins
#' @param relmarg Number insicating the fraction of the plot dedicated to the reliability score
#' @param classcale Number indicating the scaling factor for the class bars
#' @param clusterOrder Character string indicating the way to order the clusters, either score, hierarchical or none
#' @param method Character string indicating the method to show cluster reliability, either reliability or silhouette
#' @param hlim Optional vector defining the limits of the horizontal reliability plot
#' @param vlim Optional vector defining the limits of the vertical reliability plot
#' @param summary Logical whether to inclide the integrated score per cluster
#' @param silent Logical, if false, the order of the sorted samples is produced as output
#' @param cex Number indicating the relative size for the summary labels
#' @param font Number indicating the font style for the summary labels
#' @param digits Integer indicating the number of digits for the summary labels
#' @param clusterSep Optional color indicating the cluster separation lines, or vector of two colors, the first one for the bars and the second one for the heatmap
#' @param clusterSepLwd Number indicating the width of the cluster separation line, or vector of two numbers, the first one for hte bar and the second one for the heatmap
#' @param ... Additional parameters to pass to the plothm function
#' @return Nothing, a plot is generated as side effect
#' @export
plotcluster <- function(dset=NULL, hdist=NULL, vdist=NULL, hcluster=NULL, vcluster=NULL, hclass=NULL, vclass=NULL, hscore=NULL, vscore=NULL, hmargins=c(1, 1, .02, .02), relmarg=.2, classcale=2, clusterOrder=c("score", "hierarchical", "none"), method=c("reliability", "silhouette"), hlim=NULL, vlim=NULL, summary=FALSE, silent=TRUE, cex=1, font=2, digits=2, clusterSep=NULL, clusterSepLwd=1, ...) {
    clusterOrder <- match.arg(clusterOrder)
    method <- match.arg(method)
    if (!is.null(clusterSep)) if (length(clusterSep)==1) clusterSep <- rep(clusterSep, 2)
    if (length(clusterSepLwd)==1) clusterSepLwd <- rep(clusterSepLwd, 2)
    if (is.null(dset) & ((is.null(hdist) | is.null(hcluster)) & (is.null(vdist) | is.null(vcluster)))) stop("Not enough data to generate plot", call.=FALSE)
    # Only heatmap
    if ((is.null(hdist) | is.null(hcluster)) & (is.null(vdist) | is.null(vcluster))) {
        if (length(hmargins) != 4) stop("Heatmap margins should be declared as a 4 elements vector", call.=FALSE)
        warning("Only a heatmap can be generated when no distance information is provided", call.=FALSE)
        par(mai=hmargins)
        plothm(dset, ...)
        if (!is.null(colnames(dset))) axis(1, 1:ncol(dset), colnames(dset), tick=FALSE, las=2, line=-.5)
        if (!is.null(rownames(dset))) axis(2, nrow(dset):1, rownames(dset), tick=FALSE, las=2, line=-.5)
        return(TRUE)
    }
    if (!is.null(hclass)) if (is.null(nrow(hclass))) hclass <- matrix(hclass, 1, length(hclass))
    if (!is.null(vclass)) if (is.null(ncol(vclass))) vclass <- matrix(vclass, length(vclass), 1)
    # ALL
    if (all(!sapply(list(dset, hdist, vdist, hcluster, vcluster), is.null))) {
        switch(method,
               "reliability"={
                   relh <- clusterReliability(hcluster, 1/(as.matrix(hdist)+1))
                   relv <- clusterReliability(vcluster, 1/(as.matrix(vdist)+1))
               },
               "silhouette"={
                   relh <- silhouette(hcluster, dist=hdist)[, 3]
                   names(relh) <- names(hcluster)
                   relv <- silhouette(vcluster, dist=vdist)[, 3]
                   names(relv) <- names(vcluster)
                   hlim <- c(min(relh), 1)
                   vlim <- c(min(relv), 1)
               })
        if (length(hscore) == length(relh)) relh <- hscore
        if (length(vscore) == length(relv)) relv <- vscore
        posh <- sortCluster(relh, hcluster, hdist, clusterOrder)
        posv <- sortCluster(relv, vcluster, vdist, clusterOrder)
        seph <- cumsum(table(factor(hcluster[posh], levels=unique(hcluster[posh]))))[-length(unique(hcluster))]+.5 # Compute the positions for cluster separations
        sepv <- cumsum(table(factor(vcluster[rev(posv)], levels=unique(vcluster[rev(posv)]))))[-length(unique(vcluster))]+.5
        opt <- (!is.null(hclass)) + (!is.null(vclass))*10 + (relmarg>0)*100
        if (opt>0) switch(as.character(opt),
                          "1"=layout(matrix(1:2, 2, 1), heights=c(nrow(hclass)*classcale/nrow(dset), 1)),
                          "10"=layout(matrix(1:2, 1, 2), widths=c(1, ncol(vclass)*classcale/ncol(dset))),
                          "100"=layout(matrix(1:4, 2, 2), widths=c(1, relmarg), heights=c(relmarg, 1)),
                          "11"=layout(matrix(1:4, 2, 2), heights=c(nrow(hclass)*classcale/nrow(dset), 1), widths=c(1, ncol(vclass)*classcale/nrow(dset))),
                          "101"=layout(matrix(1:6, 3, 2), widths=c(1, relmarg), heights=c(relmarg, nrow(hclass)*classcale/nrow(dset), 1)),
                          "110"=layout(matrix(1:66, 2, 3), widths=c(1, ncol(vclass)*classcale/ncol(dset), relmarg), heights=c(relmarg, 1)),
                          "111"=layout(matrix(1:9, 3, 3), heights=c(relmarg, nrow(hclass)*classcale/nrow(dset), 1), widths=c(1, ncol(vclass)*classcale/ncol(dset), relmarg))
        )
        if (relmarg>0) {
            par(mai=c(hmargins[3], hmargins[2], .1, hmargins[4]))
            col <- rainbow(length(unique(hcluster)), .4, .9, start=0, end=.8)[match(hcluster[posh], unique(hcluster))]
            if (is.null(hlim)) hlim <- range(relh)
            plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", ylim=hlim, xlim=c(.5, length(relh)+.5), xaxs="i")
            for(i in 1:length(relh)) rect(i-.5, 0, i+.5, relh[posh][i], col=col[i], border=NA)
            if (summary) {
                c1 <- factor(col, levels=unique(col))
                tmp1 <- split(relh[posh], c1)
                tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=hlim)
                tmp3 <- round(tmp2, digits)
                text(cumsum(table(c1))-table(c1)/2, hlim[2]-diff(hlim)/2, tmp3, cex=cex, font=font)
            }
        }
        if (!is.null(hclass)) {
            par(mai=hmargins[c(3, 2, 3, 4)])
            image(1:ncol(hclass), 1:nrow(hclass), matrix(1:length(hclass), ncol(hclass), nrow(hclass)), col=filterColMatrix(filterRowMatrix(t(hclass), posh), nrow(hclass):1), axes=FALSE, xlab="", ylab="")
            if (!is.null(clusterSep)) abline(v=seph, col=clusterSep[1], lwd=clusterSepLwd[1])
            box()
        }
        par(mai=hmargins)
        plothm(filterColMatrix(filterRowMatrix(dset, posv), posh), ...)
        if (!is.null(clusterSep)) abline(v=seph, h=sepv, col=clusterSep[2], lwd=clusterSepLwd[2])
        if (!is.null(colnames(dset))) axis(1, 1:ncol(dset), colnames(dset)[posh], tick=FALSE, las=2, line=-.5)
        if (!is.null(rownames(dset))) axis(2, nrow(dset):1, rownames(dset)[posv], tick=FALSE, las=2, line=-.5)
        nn <- (!is.null(hclass)) + (relmarg>0)
        if (!is.null(vclass)) {
            if (nn>0) {
                for (i in 1:nn) {
                    par(mai=c(0, 0, 0, 0))
                    plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")
                }
            }
            par(mai=hmargins[c(1, 4, 3, 4)])
            image(1:ncol(vclass), 1:nrow(vclass), t(matrix(1:length(vclass), nrow(vclass), ncol(vclass))), col=filterColMatrix(filterRowMatrix(vclass, posv), ncol(vclass):1), axes=FALSE, xlab="", ylab="")
            if (!is.null(clusterSep)) abline(h=sepv, col=clusterSep[1], lwd=clusterSepLwd[1])
            box()
        }
        if (relmarg>0) {
            if (nn>0) {
                for (i in 1:nn) {
                    par(mai=c(0, 0, 0, 0))
                    plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")
                }
            }
            par(mai=c(hmargins[1], hmargins[4], hmargins[3], .1))
            col <- rainbow(length(unique(vcluster)), .4, .9, start=0, end=.8)[match(vcluster[posv], unique(vcluster))]
            if (is.null(vlim)) vlim <- range(relv)
            plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=vlim, ylim=c(.5, length(relv)+.5), yaxs="i")
            for(i in 1:length(relv)) rect(0, i-.5, relv[rev(posv)][i], i+.5, col=rev(col)[i], border=NA)
            if (summary) {
                c1 <- factor(rev(col), levels=unique(rev(col)))
                tmp1 <- split(relv[rev(posv)], c1)
                tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=vlim)
                tmp3 <- round(tmp2, digits)
                text(vlim[2]-diff(vlim)/2, cumsum(table(c1))-table(c1)/2, tmp3, cex=cex, font=font)
            }

        }
        if (silent) return(TRUE)
        return(list(x=posh, y=posv))
    }
    # Horizontal cluster and heatmap
    if (all(!sapply(list(dset, hdist, hcluster), is.null))) {
        switch(method,
               "reliability"={relh <- clusterReliability(hcluster, 1/(as.matrix(hdist)+1))},
               "silhouette"={
                   relh <- silhouette(hcluster, dist=hdist)[, 3]
                   names(relh) <- names(hcluster)
                   hlim <- c(min(relh), 1)
               })
        if (length(hscore) == length(relh)) relh <- hscore
        posh <- sortCluster(relh, hcluster, hdist, clusterOrder)
        seph <- cumsum(table(factor(hcluster[posh], levels=unique(hcluster[posh]))))[-length(unique(hcluster))]+.5 # Compute the positions for cluster separations
        if (is.null(hclass)) {
            if (relmarg>0) {
                layout(matrix(1:2, 2, 1), heights=c(relmarg, 1))
                par(mai=c(hmargins[3], hmargins[2], .1, hmargins[4]))
                col <- rainbow(length(unique(hcluster)), .4, .9, start=0, end=.8)[match(hcluster[posh], unique(hcluster))]
                if (is.null(hlim)) hlim <- range(relh)
                plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", ylim=hlim, xlim=c(.5, length(relh)+.5), xaxs="i")
                for(i in 1:length(relh)) rect(i-.5, 0, i+.5, relh[posh][i], col=col[i], border=NA)
                if (summary) {
                    c1 <- factor(col, levels=unique(col))
                    tmp1 <- split(relh[posh], c1)
                    tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=hlim)
                    tmp3 <- round(tmp2, digits)
                    text(cumsum(table(c1))-table(c1)/2, hlim[2]-diff(hlim)/2, tmp3, cex=cex, font=font)
                }
            }
            par(mai=hmargins)
            plothm(filterColMatrix(dset, posh), ...)
            if (!is.null(clusterSep)) abline(v=seph, col=clusterSep[2], lwd=clusterSepLwd[2])
            if (!is.null(colnames(dset))) axis(1, 1:ncol(dset), colnames(dset)[posh], tick=FALSE, las=2, line=-.5)
            if (!is.null(rownames(dset))) axis(2, nrow(dset):1, rownames(dset)[posv], tick=FALSE, las=2, line=-.5)
        }
        else {
            if (relmarg>0) {
                layout(matrix(1:3, 3, 1), heights=c(relmarg, nrow(hclass)*classcale/nrow(dset), 1))
                par(mai=c(hmargins[3], hmargins[2], .1, hmargins[4]))
                col <- rainbow(length(unique(hcluster)), .4, .9, start=0, end=.8)[match(hcluster[posh], unique(hcluster))]
                if (is.null(hlim)) hlim <- range(relh)
                plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", ylim=hlim, xlim=c(.5, length(relh)+.5), xaxs="i")
                for(i in 1:length(relh)) rect(i-.5, 0, i+.5, relh[posh][i], col=col[i], border=NA)
                if (summary) {
                    c1 <- factor(col, levels=unique(col))
                    tmp1 <- split(relh[posh], c1)
                    tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=hlim)
                    tmp3 <- round(tmp2, digits)
                    text(cumsum(table(c1))-table(c1)/2, hlim[2]-diff(hlim)/2, tmp3, cex=cex, font=font)
                }
            }
            else {
                layout(matrix(1:2, 2, 1), heights=c(nrow(hclass)*classcale/nrow(dset), 1))
            }
            par(mai=hmargins[c(3, 2, 3, 4)])
            image(1:ncol(hclass), 1:nrow(hclass), matrix(1:length(hclass), ncol(hclass), nrow(hclass)), col=filterColMatrix(filterRowMatrix(t(hclass), posh), nrow(hclass):1), axes=FALSE, xlab="", ylab="")
            if (!is.null(clusterSep)) abline(v=seph, col=clusterSep[1], lwd=clusterSepLwd[1])
            box()
            par(mai=hmargins)
            plothm(filterColMatrix(dset, posh), ...)
            if (!is.null(clusterSep)) abline(v=seph, col=clusterSep[2], lwd=clusterSepLwd[2])
            if (!is.null(colnames(dset))) axis(1, 1:ncol(dset), colnames(dset)[posh], tick=FALSE, las=2, line=-.5)
            if (!is.null(rownames(dset))) axis(2, nrow(dset):1, rownames(dset)[posv], tick=FALSE, las=2, line=-.5)
        }
        if (silent) return(TRUE)
        return(list(x=posh))
    }
    # Vertical cluster and heatmap
    if (all(!sapply(list(dset, vdist, vcluster), is.null))) {
        switch(method,
               "reliability"={relv <- clusterReliability(vcluster, 1/(as.matrix(vdist)+1))},
               "silhouette"={
                   relv <- silhouette(vcluster, dist=vdist)[, 3]
                   names(relv) <- names(vcluster)
                   vlim <- c(min(relv), 1)
               })
        if (length(vscore) == length(relv)) relv <- vscore
        posv <- sortCluster(relv, vcluster, vdist, clusterOrder)
        sepv <- cumsum(table(factor(vcluster[rev(posv)], levels=unique(vcluster[rev(posv)]))))[-length(unique(vcluster))]+.5
        if (is.null(vclass)) {
            if (relmarg>0) {
                layout(matrix(2:1, 1, 2), widths=c(1, relmarg))
                par(mai=c(hmargins[1], hmargins[4], hmargins[3], .1))
                col <- rainbow(length(unique(vcluster)), .4, .9, start=0, end=.8)[match(vcluster[posv], unique(vcluster))]
                if (is.null(vlim)) vlim <- range(relv)
                plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=vlim, ylim=c(.5, length(relv)+.5), yaxs="i")
                for(i in 1:length(relv)) rect(0, i-.5, relv[rev(posv)][i], i+.5, col=rev(col)[i], border=NA)
                if (summary) {
                    c1 <- factor(rev(col), levels=unique(rev(col)))
                    tmp1 <- split(relv[rev(posv)], c1)
                    tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=vlim)
                    tmp3 <- round(tmp2, digits)
                    text(vlim[2]-diff(vlim)/2, cumsum(table(c1))-table(c1)/2, tmp3, cex=cex, font=font)
                }
            }
            par(mai=hmargins)
            plothm(filterRowMatrix(dset, posv), ...)
            if (!is.null(clusterSep)) abline(h=sepv, col=clusterSep[2], lwd=clusterSepLwd[2])
            if (!is.null(colnames(dset))) axis(1, 1:ncol(dset), colnames(dset)[posh], tick=FALSE, las=2, line=-.5)
            if (!is.null(rownames(dset))) axis(2, nrow(dset):1, rownames(dset)[posv], tick=FALSE, las=2, line=-.5)
        }
        else {
            if (relmarg>0) {
                layout(matrix(3:1, 1, 3), widths=c(1, ncol(vclass)*classcale/nrow(dset), relmarg))
                par(mai=c(hmargins[1], hmargins[4], hmargins[3],.1))
                col <- rainbow(length(unique(vcluster)), .4, .9, start=0, end=.8)[match(vcluster[posv], unique(vcluster))]
                if (is.null(vlim)) vlim <- range(relv)
                plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=vlim, ylim=c(.5, length(relv)+.5), yaxs="i")
                for(i in 1:length(relv)) rect(0, i-.5, relv[rev(posv)][i], i+.5, col=rev(col)[i], border=NA)
                if (summary) {
                    c1 <- factor(rev(col), levels=unique(rev(col)))
                    tmp1 <- split(relv[rev(posv)], c1)
                    tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=vlim)
                    tmp3 <- round(tmp2, digits)
                    text(vlim[2]-diff(vlim)/2, cumsum(table(c1))-table(c1)/2, tmp3, cex=cex, font=font)
                }
            }
            else {
                layout(matrix(2:1, 1, 2), widths=c(1, nrow(hclass)*classcale/nrow(dset)))
            }
            par(mai=hmargins[c(1, 4, 3, 4)])
            image(1:ncol(vclass), 1:nrow(vclass), t(matrix(1:length(vclass), nrow(vclass), ncol(vclass))), col=filterColMatrix(filterRowMatrix(vclass, rev(posv)), ncol(vclass):1), axes=FALSE, xlab="", ylab="")
            if (!is.null(clusterSep)) abline(h=sepv, col=clusterSep[1], lwd=clusterSepLwd[1])
            box()
            par(mai=hmargins)
            plothm(filterRowMatrix(dset, posv), ...)
            if (!is.null(clusterSep)) abline(h=sepv, col=clusterSep[2], lwd=clusterSepLwd[2])
            if (!is.null(colnames(dset))) axis(1, 1:ncol(dset), colnames(dset)[posh], tick=FALSE, las=2, line=-.5)
            if (!is.null(rownames(dset))) axis(2, nrow(dset):1, rownames(dset)[posv], tick=FALSE, las=2, line=-.5)
        }
        if (silent) return(TRUE)
        return(list(y=posv))
    }
    # Horizontal only
    if (all(!sapply(list(hdist, hcluster), is.null))) {
        switch(method,
               "reliability"={relh <- clusterReliability(hcluster, 1/(as.matrix(hdist)+1))},
               "silhouette"={
                   relh <- silhouette(hcluster, dist=hdist)[, 3]
                   names(relh) <- names(hcluster)
                   hlim <- c(min(relh), 1)
               })
        if (length(hscore) == length(relh)) relh <- hscore
        posh <- sortCluster(relh, hcluster, hdist, clusterOrder)
        seph <- cumsum(table(factor(hcluster[posh], levels=unique(hcluster[posh]))))[-length(unique(hcluster))]+.5 # Compute the positions for cluster separations
        if (is.null(hclass)) {
            par(mai=c(hmargins[3], hmargins[4], .1, hmargins[4]))
            col <- rainbow(length(unique(hcluster)), .4, .9, start=0, end=.8)[match(hcluster[posh], unique(hcluster))]
            if (is.null(hlim)) hlim <- range(relh)
            plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", ylim=hlim, xlim=c(.5, length(relh)+.5), xaxs="i")
            for(i in 1:length(relh)) rect(i-.5, 0, i+.5, relh[posh][i], col=col[i], border=NA)
            if (summary) {
                c1 <- factor(col, levels=unique(col))
                tmp1 <- split(relh[posh], c1)
                tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=hlim)
                tmp3 <- round(tmp2, digits)
                text(cumsum(table(c1))-table(c1)/2, hlim[2]-diff(hlim)/2, tmp3, cex=cex, font=font)
            }
        }
        else {
            layout(matrix(1:2, 2, 1), heights=c(1, nrow(hclass)/8))
            par(mai=c(hmargins[3], hmargins[4], .1, hmargins[4]))
            col <- rainbow(length(unique(hcluster)), .4, .9, start=0, end=.8)[match(hcluster[posh], unique(hcluster))]
            if (is.null(hlim)) hlim <- range(relh)
            plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", ylim=hlim, xlim=c(.5, length(relh)+.5), xaxs="i")
            for(i in 1:length(relh)) rect(i-.5, 0, i+.5, relh[posh][i], col=col[i], border=NA)
            if (summary) {
                c1 <- factor(col, levels=unique(col))
                tmp1 <- split(relh[posh], c1)
                tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=hlim)
                tmp3 <- round(tmp2, digits)
                text(cumsum(table(c1))-table(c1)/2, hlim[2]-diff(hlim)/2, tmp3, cex=cex, font=font)
            }
            par(mai=hmargins[c(3, 4, 3, 4)])
            image(1:ncol(hclass), 1:nrow(hclass), matrix(1:length(hclass), ncol(hclass), nrow(hclass)), col=filterColMatrix(filterRowMatrix(t(hclass), posh), nrow(hclass):1), axes=FALSE, xlab="", ylab="")
            if (!is.null(clusterSep)) abline(v=seph, col=clusterSep[1], lwd=clusterSepLwd[1])
            box()
        }
        if (silent) return(TRUE)
        return(list(x=posh))
    }
    # Vertical only
    switch(method,
           "reliability"={relv <- clusterReliability(vcluster, 1/(as.matrix(vdist)+1))},
           "silhouette"={
               relv <- silhouette(vcluster, dist=vdist)[, 3]
               names(relv) <- names(vcluster)
               vlim <- c(min(relv), 1)
           })
    if (length(vscore) == length(relv)) relv <- vscore
    posv <- sortCluster(relv, vcluster, vdist, clusterOrder)
    sepv <- cumsum(table(factor(vcluster[rev(posv)], levels=unique(vcluster[rev(posv)]))))[-length(unique(vcluster))]+.5
    if (is.null(vclass)) {
        par(mai=c(hmargins[3], hmargins[4], .1, hmargins[4]))
        col <- rainbow(length(unique(vcluster)), .4, .9, start=0, end=.8)[match(vcluster[posv], unique(vcluster))]
        if (is.null(vlim)) vlim <- range(relv)
        plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", ylim=vlim, xlim=c(.5, length(relv)+.5), xaxs="i")
        for(i in 1:length(relv)) rect(i-.5, 0, i+.5, relv[rev(posv)][i], col=rev(col)[i], border=NA)
        if (summary) {
            c1 <- factor(rev(col), levels=unique(rev(col)))
            tmp1 <- split(relv[rev(posv)], c1)
            tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=vlim)
            tmp3 <- round(tmp2, digits)
            text(vlim[2]-diff(vlim)/2, cumsum(table(c1))-table(c1)/2, tmp3, cex=cex, font=font)
        }
    }
    else {
        hclass <- t(vclass)
        layout(matrix(1:2, 2, 1), heights=c(1, nrow(hclass)/8))
        par(mai=c(hmargins[3], hmargins[4], .1, hmargins[4]))
        col <- rainbow(length(unique(vcluster)), .4, .9, start=0, end=.8)[match(vcluster[posv], unique(vcluster))]
        if (is.null(vlim)) vlim <- range(relv)
        plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", ylim=vlim, xlim=c(.5, length(relv)+.5), xaxs="i")
        for(i in 1:length(relv)) rect(i-.5, 0, i+.5, relv[rev(posv)][i], col=rev(col)[i], border=NA)
        if (summary) {
            c1 <- factor(rev(col), levels=unique(rev(col)))
            tmp1 <- split(relv[rev(posv)], c1)
            tmp2 <- sapply(tmp1, function(x, xlim) 1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim), xlim=vlim)
            tmp3 <- round(tmp2, digits)
            text(vlim[2]-diff(vlim)/2, cumsum(table(c1))-table(c1)/2, tmp3, cex=cex, font=font)
        }
        par(mai=hmargins[c(3, 4, 3, 4)])
        image(1:ncol(hclass), 1:nrow(hclass), matrix(1:length(hclass), ncol(hclass), nrow(hclass)), col=filterColMatrix(filterRowMatrix(t(hclass), posv), nrow(hclass):1), axes=FALSE, xlab="", ylab="")
        if (!is.null(clusterSep)) abline(h=sepv, col=clusterSep[1], lwd=clusterSepLwd[1])
        box()
    }
    if (silent) return(TRUE)
    return(list(y=posv))
}

#' Sort the samples after cluster analysis
#'
#' This function generate an index to sort the samples according a cluster structure and distance
#' @param rel Reliability score as computed by \code{clusterReliability}
#' @param cluster Vector of cluster membership
#' @param dis Distance object
#' @param clusterOrder Character string indicating the method for sorting the clusters, either score or hierarchical
#' @return Numeric vector of interger positions
sortCluster <- function(rel, cluster, dis, clusterOrder=c("score", "hierarchical", "none")) {
    cl <- cluster
    if (clusterOrder != "none") {
        dis <- as.matrix(dis)
        cdis <- sapply(unique(cluster), function(i, cluster, dis, rel) {
            sapply(unique(cluster), function(ii, i, cluster, dis, rel) {
                ii <- which(cluster==ii)
                i <- which(cluster==i)
                tmp <- dis[ii, ][, i]
                wt <- matrix(rep(rel[ii], length(i))*rep(rel[i], each=length(ii)), length(ii), length(i))
                sum(tmp*wt/sum(wt))
            }, i=i, cluster=cluster, dis=dis, rel=rel)
        }, cluster=cluster, dis=dis, rel=rel)
        colnames(cdis) <- rownames(cdis) <- unique(cluster)
        pos <- hclust(as.dist(cdis))$order
        cl <- unique(cluster)[match(match(cluster, unique(cluster)), pos)]
        if (clusterOrder=="hierarchical") {
            pos <- split(1:length(rel), cluster)
            tmp <- lapply(pos, function(i, dis) {
                match(1:length(i), hclust(as.dist(dis[i, ][, i]))$order)
            }, dis=dis)
            rel[unlist(pos, use.names=FALSE)] <- unlist(tmp, use.names=FALSE)
        }
    }
    order(cl, -rel)
}

#' Optimal cluster
#'
#' This function generate an optimal cluster structure based on cluster reliability score analysis
#'
#' @param dis Distance object
#' @param step Integer indicating the incremental number of clusters to add in each iteration
#' @param cores Maximum number of CPU cores to use
#' @return Vector of cluster membership
#' @export
optimalCluster <- function(dis, step=5, cores=1) {
    nmax <- round(attr(dis, "Size")/3)
    if (nmax<(step*2)) {
        clus <- mclapply(2:(2*step), function(i, dis) pam(dis, i, diss=TRUE, cluster.only=TRUE), dis=dis, mc.cores=min(cores, step*2-1))
        return(clus[[which.max(unlist(clusterReliability(clus, 1/as.matrix(dis), method="global"), use.names=FALSE))]])
    }
    clus1 <- mclapply(2:step, function(i, dis) pam(dis, i, diss=TRUE, cluster.only=TRUE), dis=dis, mc.cores=min(cores, step-1))
    pointer <- step
    clus2 <- mclapply((pointer+1):(pointer+step), function(i, dis) pam(dis, i, diss=TRUE, cluster.only=TRUE), dis=dis, mc.cores=min(cores, step))
    cr <- clusterReliability(c(clus1, clus2), 1/as.matrix(dis), method="global")
    while(which.max(cr)>(length(cr)-step) & pointer < nmax) {
        clus1 <- c(clus1, clus2)
        pointer <- pointer + step
        clus2 <- mclapply((pointer+1):(pointer+step), function(i, dis) pam(dis, i, diss=TRUE, cluster.only=TRUE), dis=dis, mc.cores=min(cores, step))
        cr <- clusterReliability(c(clus1, clus2), 1/as.matrix(dis), method="global")
    }
    clus1[[which.max(cr)]]
}
