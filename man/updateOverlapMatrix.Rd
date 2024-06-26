% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processNewRamp.R
\name{updateOverlapMatrix}
\alias{updateOverlapMatrix}
\title{Update the overlap matrix for store in the shiny app directory}
\usage{
updateOverlapMatrix(min_analyte, overlapmethod, together)
}
\arguments{
\item{min_analyte}{a int that specifies the minimum of analytes the
pathway should have to be considered compute for overlap matrix}

\item{overlapmethod}{a string that specifies the way to compute overlap matrix,
must be 'balanced' or 'weighted'}

\item{together}{a boolean value to compute overlap matrix for
gene/metabolites separatly or together}
}
\description{
Update the overlap matrix for store in the shiny app directory
}
