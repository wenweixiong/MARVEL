% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Script_PLATE_00_rMATS_00_CreateFeatureData.R
\name{Preprocess_rMATS}
\alias{Preprocess_rMATS}
\title{Preprocess rMATS files for MARVEL}
\usage{
Preprocess_rMATS(file, GTF, EventType)
}
\arguments{
\item{file}{Data frame. Tab-delimited file output from rMATS. Typically the file is named \code{fromGTF.*.txt}.}

\item{GTF}{Data frame. The same GTF file used in the rMATS step.}

\item{event_type.}{Indicate the type of splicing event. Options are \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"}.}
}
\value{
An object of class data frame that may be used as input for MARVEL
}
\description{
Preprocess rMATS coordinates for MARVEL input and annotate the gene type for each gene using the GTF provided.
}
