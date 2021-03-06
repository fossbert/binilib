% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{find_var_genes}
\alias{find_var_genes}
\title{Finding the most variable features of a matrix}
\usage{
find_var_genes(expmat, n = 2000, method = c("mad", "sd"), center = FALSE)
}
\arguments{
\item{expmat}{numerical matrix, typically a gene expression matrix}

\item{n}{number of most variable features to be returned, default 2000}

\item{method}{method to determine 'variable', defaults to robust measure (mad)}

\item{center}{whether to median center the resulting matrix, defaults to FALSE}
}
\value{
new matrix with most variable genes
}
\description{
Finding the most variable features of a matrix
}
