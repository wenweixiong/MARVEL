% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/Script_PLATE_00_rMATS_05_CreateFeatureData_A3SS.R
\name{Preprocess_rMATS.A3SS}
\alias{Preprocess_rMATS.A3SS}
\title{Preprocess rMATS alternative 3' splice site (A3SS) coordinates}
\usage{
Preprocess_rMATS.A3SS(file, GTF)
}
\arguments{
\item{file}{Data frame. Tab-delimited file output from rMATS. Typically the file is named \code{fromGTF.A3SS.txt}.}

\item{GTF}{Data frame. The same GTF file used in the rMATS step.}
}
\value{
An object of class data frame that may be used as input for MARVEL
}
\description{
Preprocess rMATS alternative 5' splice site (A3SS) coordinates for MARVEL input and annotate the gene type for each gene using the GTF provided.
}
