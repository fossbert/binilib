#' Information on gene lengths.
#'
#' As derived from TxDb.Hsapiens.UCSC.hg19.knownGene and TxDb.Mmusculus.UCSC.mm9.knownGene annotation
#'
#' @format List with two elements containing human and mouse data, respectively, each with:
#' \describe{
#' \item{Symbol}{named integer vector of gene lengths, where names represent gene symbols}
#' \item{Entrez}{named integer vector of gene lengths, where names represent Entrez IDs}
#' }
#' @source \url{http://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html}
#' @source \url{http://bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm9.knownGene.html}
"geneLengths"

#' Human and Mouse Homology data
#'
#' Information on 16469 genes for which there is exactly a 1:1 mapping from a mouse gene id to a human gene id.
#'
#' @format List with two elements containing gene \strong{symbol} and \strong{entrez id} data, respectively
#' \describe{
#' \item{Symbol}{named character vector of mouse gene symbols, where names represent human symbols}
#' \item{Entrez}{named character vector of mouse gene entrez ids, where names represent human entrez ids}
#' }
#' @source \url{http://www.informatics.jax.org/homology.shtml}
"homologyInfo"

#' Gene expression signatures used to classify PDA.
#'
#' Derived from Collisson et al., Moffitt et al. and Bailey et al.
#'
#' @format List with 11 elements containing gene signatures (symbol identifiers)
#'
"pdaClassifiers"
