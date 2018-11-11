#' Extract the top up and down regulated MR per sample
#'
#' @param vpres numeric matrix of NES as returned by the viper function
#' @param nn integer of number of top MR to consider for each regulon
#' @param direction character, which tail should be returned, defaults to both up and down
#' @param make_regulon logical, whether to set class of return object to "regulon"
#' @param reverse, logical, whether to reverse the top MR - useful for drug screening results and downstream applications
#' @param verbose logical, whether to message progress report
#' @return list object with top dysregulated MR per sample
#' @export

vpres2regulon <- function(vpres,
                          nn = 25,
                          direction = c('both', 'up', 'down'),
                          make_regulon = FALSE,
                          reverse = FALSE,
                          verbose = TRUE) {

    dir <- match.arg(direction)

    switch(dir,
           both = {idx <- c(1:nn, (nrow(vpres)-(nn-1)):nrow(vpres))
           tfm <- rep(c(1,-1), each = nn) },
           up = {idx <- 1:nn
           tfm <- rep(1, nn)},
           down = {idx <- (nrow(vpres)-(nn-1)):nrow(vpres)
           tfm <- rep(-1, nn)}
    )

    if(reverse) tfm <- tfm * -1

    if(verbose){
        message('Building regulon ...')
        pb <- txtProgressBar(min = 0, max = ncol(vpres), style = 3)
    } else pb <- NULL

    tmp <- lapply(1:ncol(vpres), function(i, vpres) {

        if(!is.null(pb)) setTxtProgressBar(pb, i)

        smp <- vpres[,i]
        names(tfm) <- rownames(vpres[order(smp, decreasing = TRUE), ])[idx]

        list(tfmode = tfm, likelihood = unname(abs(tfm)))

    }, vpres = vpres)
    names(tmp) <- colnames(vpres)

    if(make_regulon) class(tmp) <- "regulon"
    return(tmp)
}

#' transform a regulon object to a data frame / tibble
#'
#' @param reg regulon object - see output of viper::aracne2regulon
#' @param annotate logical, whether to add a column on target type in ouptu
#' @return a tibble with targets and metrics (MOR/likelihood)
#' @export

reg2tibble <- function(reg, annotate = TRUE) {

        tmp <- tibble::tibble(source = rep(names(reg), purrr::map_int(reg, ~ length(.$tfmode))),
               target = unlist(purrr::map(reg, ~ names(.$tfmode)), use.names = FALSE),
               mor = unlist(purrr::map(reg, ~ .$tfmode), use.names=FALSE),
               likelihood = unlist(purrr::map(reg, ~ .$likelihood), use.names=FALSE)
               ) %>%
            dplyr::arrange(source, target, mor)

        if(annotate) {
            data("regulatorNames")

            tmp2 <- tibble::tibble(target = unlist(regulators, use.names = FALSE),
                                   type = rep(names(regulators), purrr::map_int(regulators, length)))

            dplyr::left_join(tmp, tmp2, by = "target") %>%
                dplyr::mutate(type = stringr::str_replace_na(type, "None"))

        } else tmp

}

#' viperSignature and nullmodel using limma
#'
#' @param eset normalized gene expression matrix with samples of interest
#' @param ref normalized gene expression matrix with reference samples
#' @param per integer, number of overall permutations
#' @param keep_all logical, whether to return signature for all samples (incl. reference samples)
#' @param verbose logical, whether to report on progress or not
#' @export

vpsig_limma <- function(eset, ref, per = 1000, keep_all = FALSE, seed = 42, verbose = TRUE) {

    set.seed(seed)
    if(keep_all) eset <- cbind(eset, ref) # put out all samples ?

    # calculate signature using moderated t-statistics
    if (verbose) {
        message("\nCalculating the signature for ", ncol(eset), " samples")
        pb <- txtProgressBar(max = ncol(eset), style = 3)
    } else pb <- NULL

    signature <- vapply(1:ncol(eset), function(i, dset, refmat){

            if (!is.null(pb)) setTxtProgressBar(pb, i)
            tmp <- eBayes(lmFit(dset[,i] - refmat))
            qnorm(tmp$p.value[,1]/2, lower.tail = FALSE) * sign(tmp$t[,1])

    }, FUN.VALUE = numeric(nrow(eset)), dset = eset, refmat = ref)
    colnames(signature) <- colnames(eset)

    # calculate nullmodel
    if (ncol(ref)<6) {
        vpnull <- NULL
        warning("Not enough samples to compute null model by sample permutation,
                gene permutation will be used instead", call.=FALSE)
    } else {

        # message methods of collapsing, number of samples and date
        if (verbose) {
            message("\nCalculating nullmodel using bootstraps of the ", ncol(ref), " reference samples")
            pb <- txtProgressBar(max = ncol(ref), style = 3)
        } else pb <- NULL


        # leave-one-out vs. bootstrap set

      vpnull <- do.call(cbind, lapply(1:ncol(ref), function(i, d22, per){

            if (!is.null(pb)) setTxtProgressBar(pb, i)

            vapply(1:per, function(i, x, y) {

                pos <- sample(ncol(y), replace=TRUE)
                tmp <- eBayes(lmFit(x - y[ ,pos]))
                qnorm(tmp$p.value[,1]/2, lower.tail = FALSE) * sign(tmp$t[,1])

            }, FUN.VALUE = numeric(nrow(d22)), x = d22[ ,i], y = d22[ ,-i]) # vapply end

        }, d22 = ref, per = ceiling(per/ncol(ref)))) # lapply end
    }

    # put together output
    vpsig <- list(signature = signature, nullmodel = vpnull)
    class(vpsig) <- "viperSignature"
    return(vpsig)

}


#' Convert regular (one-tailed) gene sets to regulons
#'
#' So they can be used with the viper function
#'
#' @param geneSets list of named lists representing gene sets of gene symbols
#' @return regulon containing one-tailed gene sets as regulators
#' @export

geneSets2regulon <- function(geneSets) {

    tmp <- lapply(geneSets, function(i){

        tm <- rep(1, length(i))
        names(tm) <- i
        list(tfmode = tm, likelihood = rep(1/length(i), length(i)))

    })
    class(tmp) <- 'regulon'
    return(tmp)
}

#' Extract msviper results
#'
#' @param mrst msviper object
#' @param padjust method for adjusting p-values for multiple hypothesis testing
#' @return a tibble with genes, NES, p-values and FDR
#' @export

res_msvip <- function(mrst, padjust = c('fdr', 'bonferroni','none')) {

    method <- match.arg(padjust, choices = c('fdr', 'bonferroni','none'))

    if(is.null(mrst$es$nes.se)){
        res <- tibble::tibble(gene = names(mrst$es$nes),
                              nes = mrst$es$nes,
                              pval = mrst$es$p.value,
                              padj = p.adjust(mrst$es$p.value, method = method))
    } else{
        res <- tibble::tibble(gene = names(mrst$es$nes),
                              nes = mrst$es$nes,
                              nes_se = mrst$es$nes.se,
                              pval = mrst$es$p.value,
                              padj = p.adjust(mrst$es$p.value, method = method))
    }

    res <- dplyr::arrange(res, desc(nes))
    return(res)

}

