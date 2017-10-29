#' Extract rows or columns from a matrix in order of interest
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
