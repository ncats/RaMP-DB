% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processNewRamp.R
\name{findAnalyteHasPathway}
\alias{findAnalyteHasPathway}
\title{Find table of analyte has pathway from given pathway IDs
Aggregate ramp Id to ramp pathway Id
GC is C or G}
\usage{
findAnalyteHasPathway(db = RaMP(), pathwayRampId, GC = "C", n = 10)
}
\arguments{
\item{pathwayRampId}{a vector of ramp Pathway ID}

\item{GC}{the analytes type that is either "C" for compound or "G" for gene}

\item{n}{minimum analytes of which pathway to considered computing overlap}
}
\value{
A list with pathway rampID as name, a vector of analytes from this pathway as content.
}
\description{
Find table of analyte has pathway from given pathway IDs
Aggregate ramp Id to ramp pathway Id
GC is C or G
}
