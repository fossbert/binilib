#' PurIST
#'
#' This function implements single-sample classification for PDAC specimen based on their raw gene expression profiles.
#' Idealy you would normalize for gene length although it does work with microarray intensities as well, so should be fine
#' with normalized counts just accounting for different library sizes.
#'
#' PMID: 31754050
#'
#' Probability key:
#' < 0.1: strongly favors classical
#' 0.1 - 0.4: likely classical
#' 0.4 - 0.5: lean classical
#' 0.5 - 0.6: lean basal-like
#' 0.6 - 0.9: likely basal-like
#' >0.9: strongly favors basal-like
#'
#' @param emat expression matrix with row (=gene) and column (=sample) names. Genes must be supplied as HUGO gene symbols
#' @param raw logical indicating on whether to return a matrix with raw expression values for all pairs for inspection (default = FALSE)
#' @return a named numeric vector of probabilities reflecting class membership or a matrix of pasted expression values
#' @export
#'
purIST <- function(emat, raw = FALSE){

    stopifnot(is.matrix(emat))
    stopifnot(!is.null(rownames(emat)))

    # Data - see Rashid et al. Supplement
    coefs <- c(-6.815, 1.994, 2.031, 1.618, 0.922, 1.059, 0.929, 2.505, 0.485)

    tsps <- list(
        c('GPR87', 'REG4'),
        c('KRT6A', 'ANXA10'),
        c('BCAR3', 'GATA6'),
        c('PTGES', 'CLDN18'),
        c('ITGA3', 'LGALS4'),
        c('C16orf74', 'DDC'),
        c('S100A2', 'SLC40A1'),
        c('KRT5', 'CLRN3')
    )

    predictors <- unlist(tsps, use.names = FALSE)

    # check for avalaibility of all predictors
    if(!all(predictors %in% rownames(emat))) stop('Could not find all genes needed for PurIST evaluation')

    if (raw) {

        tmp <- apply(round(emat, 3), 2, function(i){

            vapply(tsps, function(j) {

                paste(i[j[1]], i[j[2]], sep = '_')

            }, FUN.VALUE = character(1))
        })

        rownames(tmp) <- vapply(tsps, paste, collapse = '_', FUN.VALUE = character(1))

        return(tmp)
    } else {

        # Create indicator matrix
        tmp <- apply(emat, 2, function(i){

            c(1, vapply(tsps, function(j) {

                as.numeric((i[j[1]] > i[j[2]])-.5 > 0)

            }, FUN.VALUE = numeric(1)))
        })
        # calculate and return probabilities
        logits <- t(tmp) %*% coefs
        probos <- exp(logits)/(1 + exp(logits))
        return(probos[,1])

    }

}
