% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting.R
\name{legend_col}
\alias{legend_col}
\title{Add a custom color gradient to a base R plot}
\usage{
legend_col(col, level, side = 4)
}
\arguments{
\item{col}{a chracter vector as returned from colorRampPalette}

\item{level}{a factor generated by cutting a continuous variable into many pieces (e.g. 100)}

\item{side}{where to draw the color gradient, default is 4, i.e. right}
}
\value{
added a color gradient to existing plot
}
\description{
Add a custom color gradient to a base R plot
}
