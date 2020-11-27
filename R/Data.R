#' Information on gene lengths.
#' As derived from TxDb.Hsapiens.UCSC.hg19.knownGene and TxDb.Mmusculus.UCSC.mm9.knownGene annotation
#'
#' @name geneLengths
#' @docType data
#' @author HCM
#' @keywords data
#' @format List with two elements containing human and mouse data, respectively, each with:
#' \describe{
#' \item{Symbol}{named integer vector of gene lengths, where names represent gene symbols}
#' \item{Entrez}{named integer vector of gene lengths, where names represent Entrez IDs}
#' }
#' @source \url{http://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html}
#' @source \url{http://bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm9.knownGene.html}
NULL

#' Human and Mouse Homology data
#' Information on 16469 genes for which there is exactly a 1:1 mapping from a mouse gene id to a human gene id.
#'
#' @name homologyInfo
#' @docType data
#' @author HCM
#' @keywords data
#' @format List with two elements containing gene \strong{symbol} and \strong{entrez id} data, respectively
#' \describe{
#' \item{Symbol}{named character vector of mouse gene symbols, where names represent human symbols}
#' \item{Entrez}{named character vector of mouse gene entrez ids, where names represent human entrez ids}
#' }
#' @source \url{http://www.informatics.jax.org/homology.shtml}
NULL

#' Information on 6068 regulatory genes used in the Califanoverse
#' See Alvarez, Nature Genetics, 2016 for details
#'
#' @name regulatorNames
#' @docType data
#' @author HCM
#' @keywords data
#' @format List with three elements containing gene symbols for the following categories
#' \describe{
#' \item{TF}{character vector transcription factor gene symbols, n = 1856}
#' \item{CoTF}{character vector of co-transcriptional regulator gene symbols, n = 672}
#' \item{Signaling}{character vector of signaling proteins, n = 3540}
#' }
#'
NULL

#' Full homology data as provided by Jax
#' Information on a more complete set of 18548 genes for which the Jax Informatics Portal
#' provides mappings between human and mouse genes. This includes human genes for which there are
#' multiple corresponding mouse genes, and vice versa, mouse genes with multiple corresponding human
#' genes.
#'
#' @name homologyTable
#' @docType data
#' @author HCM
#' @keywords data
#' @format List with two elements (tibble) containing gene \strong{symbol} and \strong{entrez id} data, respectively
#' @source \url{http://www.informatics.jax.org/homology.shtml}
NULL

