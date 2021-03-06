% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSEA_null.R
\docType{methods}
\name{gsea_null}
\alias{gsea_null}
\alias{gsea_null,gsea1-method}
\alias{gsea_null,gsea2-method}
\alias{gsea_null,}
\alias{gsea2-method}
\title{GSEA null model}
\usage{
gsea_null(gsea_obj, ...)

\S4method{gsea_null}{gsea1}(
  gsea_obj,
  w = 1,
  perm = 1000,
  seed = 42,
  analytical = FALSE,
  verbose = TRUE
)

\S4method{gsea_null}{gsea2}(gsea_obj, w = 1, seed = 42, perm = 1000, verbose = TRUE)
}
\arguments{
\item{gsea_obj}{a gsea object as returned by gsea1T or gsea2T/gsea_regulon}

\item{...}{Additional parameters added to keep compatibility}

\item{w}{numeric indicating the weights used for GSEA (defaults to 1)}

\item{perm}{integer indicating the number of permutations to carry out (defaults to 1000)}

\item{seed}{integer for random number generation during the permutation process (defaults to 42)}

\item{analytical}{logical whether to derive p-values analytically, i.e. from two-tailed test using a normal distribution (defaults to FALSE, i.e. inference according to Subramanian et al.)}
}
\value{
gsea_obj with information on nullmodel added in extra slots
}
\description{
This function generates a GSEA null distribution of enrichment scores based on gene shuffling, i.e.
the positions of genes in the signature are permuted rather than phenotype labels. It will add a
normalized enrichment score (NES) to the gsea object as described by Subramanian et al. and statistical
inference is carried out using one-tailed tests on the appropriate (positive/negative)
side of the null distribution. For two-tailed GSEA, the null model is calculated accordingly and
separate NES and p-values are provided for each gene set.
}
