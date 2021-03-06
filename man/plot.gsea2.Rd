% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_gsea.R
\name{plot.gsea2}
\alias{plot.gsea2}
\title{Plot 2-tailed enrichment results}
\usage{
\method{plot}{gsea2}(
  gsea2,
  color = c("cornflowerblue", "salmon"),
  legend = TRUE,
  legNames = c("Targets", "pos", "neg"),
  main = "",
  plotSignature = FALSE,
  signatureNames = NULL,
  signatureType = NULL,
  ledge = TRUE,
  fixY = NULL,
  ...
)
}
\arguments{
\item{gsea2}{gsea2 object as returned from \code{gsea2T}}

\item{color}{Vector of two components indicating the colors for the negative and positive parts of the regulon}

\item{legend}{Logical on whether to draw a legend explaining the two parts of the regulon}

\item{legNames}{Character, defaults to Targets, pos and neg}

\item{main}{Character vector given as argument for convenience}

\item{plotSignature}{Logical whether to include a graph on the signature (e.g. t-statistics)}

\item{signatureNames}{Character vector indicating the conditions being compared. Defaults to NULL. Otherwise a charater vector of length two.}

\item{signatureType}{Character to indicate the type of gene level statistic (e.g. t-statistic).}

\item{ledge}{logical whether to plot the leading edge line}

\item{fixY}{numeric vector specyfying y-axis limits, defaults to NULL
Defaults to NULL which will yield a 'signature score' label.}

\item{...}{Given for compatibility to the plot generic function}
}
\value{
Nothing, a plot is generated in the default output device
}
\description{
This function generates a plot of the GSEA running sums of two gene sets on a given signature.
}
