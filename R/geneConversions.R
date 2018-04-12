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
    if(!base::exists("homologyInfo", where = ".GlobalEnv")) data("homologyInfo")
    common <- intersect(mgenes, homologyInfo[[type]])
    idx <- match(common, homologyInfo[[type]])

    if (length(common) == 0) stop("Sorry, could not align gene identifiers, please double check.")

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
#' @param verbose logical whether to message progress reports.
#' @return a character vector of gene symbols, beware of NAs
#' @export

any2symbol <- function(gene_ids, verbose = TRUE) {

    gene_ids <- as.character(gene_ids) # safety net, factors and integers cause trouble !
    if(!base::exists('geneInfo', where = ".GlobalEnv")) data("geneInfo")
    idx <- Position(function(i) any(gene_ids %in% i), geneInfo)

    if(is.na(idx)) stop("Sorry, could not find any of the provided gene ids !", call. = FALSE)

    # report overlaps
    if (verbose) {
        message("\nFound ", sum(gene_ids %in% geneInfo[[idx]]), " out of ", length(gene_ids), " gene identifiers")
        message("\nConverting ", tolower(names(geneInfo)[idx]), " to gene symbols")
    }

    symbol <- geneInfo[["SYMBOL"]]
    names(symbol) <- geneInfo[[idx]]
    unname(symbol[gene_ids])

}

#' Convert Symbols or Ensembl gene ids to Entrez Ids
#'
#' Two reasons for NA return: the external gene id could not be found in the reference (org.Hs.Db)
#' or there is no Entrez Id for that particular symbol or Ensembl Id, respectively.
#'
#' @param gene_ids a character vector of gene ids from the Hugo gene sybmol or Ensembl annotation, respectively
#' @param verbose logical whether to message progress reports.
#' @return a character vector of gene symbols, beware of NAs
#' @export

any2entrez <- function(gene_ids, verbose = TRUE) {

    gene_ids <- as.character(gene_ids) # safety net, factors and integers cause trouble !
    if(!base::exists('geneInfo', where = ".GlobalEnv")) data("geneInfo")
    idx <- Position(function(i) any(gene_ids %in% i), geneInfo)

    if(is.na(idx)) stop("Sorry, could not find any of the provided gene ids !", call. = FALSE)

    # report overlaps
    if (verbose) {
        message("\nFound ", sum(gene_ids %in% geneInfo[[idx]]), " out of ", length(gene_ids), " gene identifiers")
        message("\nConverting ", tolower(names(geneInfo)[idx]), " to Entrez Ids")
    }

    entrez <- geneInfo[["ENTREZID"]]
    names(entrez) <- geneInfo[[idx]]
    unname(entrez[gene_ids])

}
