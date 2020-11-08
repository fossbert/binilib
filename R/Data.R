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

#' Gene expression signatures used to classify PDA.
#' Derived from Collisson et al., Moffitt et al., Chan-Seng-Yue et al. and Bailey et al.
#'
#' @name pdaClassifiers
#' @docType data
#' @author HCM
#' @keywords data
#'
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

#' Modules c1 through c6 and hallmark gene sets from the MSigDb v6.2
#' Downloaded from the Broad Institute on October 2018
#'
#' @name msigdb
#' @docType data
#' @author HCM
#' @keywords data
#' @format List with 12 elements containing human gene symbols
#' \describe{
#' \item{c1.all}{hallmark gene sets, n = 326}
#' \item{c2.cgp}{canonical pathways, n = 3433}
#' \item{c2.cp}{transcription factor targets, n = 1329}
#' \item{c3.mir}{micro RNA targets, n = 221}
#' \item{c3.tft}{transcription factor targets, n = 615}
#' \item{c4.cgn}{cancer gene neighborhoods, n = 427}
#' \item{c4.cm}{cancer modules factor, n = 431}
#' \item{c5.bp}{GO biological processes, n = 4436}
#' \item{c5.cc}{GO cellular components, n = 580}
#' \item{c5.mf}{GO molecular function, n = 901}
#' \item{c6.all}{oncogenic signatures, n = 189}
#' \item{h.all}{hallmark gene sets, n = 50}
#' }
NULL

#' Murine version of the Selection of MSigDb gene sets v6.2
#' Downloaded from the Broad Institute Octobre 2018 and
#' converted using information on gene homolgy from the
#' Jax Bioinformatics portal
#'
#' @name msigdbMouse
#' @docType data
#' @author HCM
#' @keywords data
#' @format List 12 elements containing mouse gene symbols
#' \describe{
#' \item{c1.all}{hallmark gene sets, n = 326}
#' \item{c2.cgp}{canonical pathways, n = 3433}
#' \item{c2.cp}{transcription factor targets, n = 1329}
#' \item{c3.mir}{micro RNA targets, n = 221}
#' \item{c3.tft}{transcription factor targets, n = 615}
#' \item{c4.cgn}{cancer gene neighborhoods, n = 427}
#' \item{c4.cm}{cancer modules factor, n = 431}
#' \item{c5.bp}{GO biological processes, n = 4436}
#' \item{c5.cc}{GO cellular components, n = 580}
#' \item{c5.mf}{GO molecular function, n = 901}
#' \item{c6.all}{oncogenic signatures, n = 189}
#' \item{h.all}{hallmark gene sets, n = 50}
#' }
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

#' Information on ligands and receptors
#' As derived from Ramilowski et al., Nature Communications, 2016. PMID: 26198319.
#'
#' @name ligandsReceptors
#' @docType data
#' @author HCM
#' @keywords data
#' @format tibble with information on ligands and their receptors
NULL


#' Information on source gene sets of MSigDb Hallmark collection
#' See Liberzon, Cell Syst, 2015 for details. PMID: 26771021
#'
#' @name hallmarkFounders
#' @docType data
#' @author HCM
#' @keywords data
#' @format list of lists: for 50 Hallmark gene sets, one list each comprising all founders
NULL


#' ARACNe regulons derived from a) 242 CUMC pancreatic specimen (PanIN/IPMN/PDA) and b) 426 ICGC
#' Canada PDA specimen (PDA primaries and metastasis)
#'
#' @name pda_regulons
#' @docType data
#' @author HCM
#' @keywords data
#' @format list of regulons, each comprising one regulon with its RP
NULL

