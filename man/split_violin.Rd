% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting.R
\name{split_violin}
\alias{split_violin}
\title{Split violin plot}
\usage{
split_violin(
  x,
  s1,
  s2,
  cols = c("cornflowerblue", "salmon"),
  ylb = "",
  rug = FALSE,
  legpos = "topleft",
  rotate_xlabels = NULL,
  ...
)
}
\arguments{
\item{x}{a numeric variable to be split in two ways}

\item{s1}{categorical variable used for the first split}

\item{s2}{categorical variable used for the second split}

\item{cols}{vector of lenght 2 indicating colors in R (hex, integer or plane name)}

\item{ylb}{default for y-axis, could be changed if you know the nature of x (e.g. Z-score)}

\item{rug}{logical, whether to plot the data points as a rug}

\item{legpos}{character, indicating the position of the legend}

\item{rotate_xlabels}{integer, if not NULL, will be used as degree rotation of x labels}
}
\value{
a plot will be send to the graphics output
}
\description{
This type of plot takes a numeric vector and splits it twice. Once using a first categorical variable with
any number of levels and then and then a second time using a second categorical variable with
exactly TWO levels. A density curve will be estimated and side-by-side violin plots will be produced.
}
