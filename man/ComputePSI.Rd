% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Script_01_0_ComputePSI.R
\name{ComputePSI}
\alias{ComputePSI}
\title{Compute Percent Spliced-in (PSI) Values}
\usage{
ComputePSI(MarvelObject, CoverageThreshold, EventType,
  IntronCountsFile = NULL, thread = NULL)
}
\arguments{
\item{MarvelObject}{S3 object generated from \code{CreateMarvelObject} function.}

\item{CoverageThreshold}{Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.}

\item{EventType}{Character string. Indicate which splicing event type to calculate the PSI values for. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.}

\item{IntronCountsFile}{Data frame containing per base coverage of introns. Only applicable when \code{EventType} set to \code{"RI"}. First column should be named \code{coord.intron} and indicate the per base intron position in the form of chr:position. Subsequent columns should contain the per base coverage for each sample. These counts can be deteted using external softwares such as Bedtools etc..}

\item{thread}{Numeric value. Set number of threads.}
}
\value{
An object of class S3 containing all the original slots as inputted by the user in addition to two new slots. \code{$SpliceFeatureValidated} contains the validated splicing event metadata. \code{$PSI} contains the computed PSI values for the validated splicing events.
}
\description{
\code{ComputePSI} computes percent spliced-in (PSI) values for a specified splicing event type.
}
\details{
This function computes the PSI values for a splicing event as specified in the \code{EventType} argument. Splicing events provided in \code{SpliceFeature} data frame will first be cross-checked against the splice junctions provided in \code{SpliceJunction} data frame. Only events whose junctions are found in \code{SpliceJunction} are retained. The formula for computing PSI is the number of junction reads supporting the included isoform divided by the total number of reads supporting both included and excluded isoforms.
}
\examples{
marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                    EventType="SE"
                    )

marvel$SpliceFeatureValidated$SE
marvel$PSI$SE[,1:5]
}
\author{
Sean Wen <sean.wenwx@gmail.com>
}
