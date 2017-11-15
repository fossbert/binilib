#' Convert mouse gene symbols or entrez ids to human counterparts
#'
#' @param eset gene expression matrix with murine expression data
#' @param type Character to indicate the gene identifier, defaults to gene symbol
#' @return Numeric matrix with human gene identifiers of the same type
#' @export

mouse2human <- function(eset, type = c('Symbol', 'Entrez')) {

    mgenes <- rownames(eset)

    type <- match.arg(type, choices = c('Symbol', 'Entrez'))

    # retrieve information on gene homologies
    data("homologyInfo")
    common <- intersect(mgenes, homologyInfo[[type]])
    idx <- match(common, homologyInfo[[type]])

    if (length(common) == 0) {
        stop("Sorry, could not align gene identifiers, please double check.")
    }

    tmp <- eset[common, ]
    rownames(tmp) <- names(homologyInfo[[type]])[idx]
    return(tmp)
}


#' Convert Entrez or Ensembl gene ids to gene symbols
#'
#' Two reasons for NA return: the external gene id could not be found in the reference (org.Hs.Db)
#' or there is no gene symbol for that particular external gene id.
#'
#' @param gene_ids a character vector of gene ids from the Ensembl or Entrez annotation, respectively
#' @return a character vector of gene symbols, beware of NAs
#' @export

any2symbol <- function(gene_ids) {

    data("geneInfo")
    idx <- Position(function(i) any(gene_ids %in% i), geneInfo)

    if(is.na(idx)) stop("Sorry, could not find any of the provided gene ids among the Ensembl and
                        Entrez reference, respectively!")

    # report overlaps
    message("\nFound ", sum(gene_ids %in% geneInfo[[idx]]), " out of ", length(gene_ids), " gene identifiers")
    message("\nConverting ", tolower(names(geneInfo)[idx]), " to gene symbols")

    symbol <- geneInfo[["SYMBOL"]]
    names(symbol) <- geneInfo[[idx]]
    unname(symbol[gene_ids])

}
