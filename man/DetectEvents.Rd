% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Script_PLATE_03_DETECT_EVENTS_0.R
\name{DetectEvents}
\alias{DetectEvents}
\title{Detect Splicing Events}
\usage{
DetectEvents(
  MarvelObject,
  min.cells = 50,
  min.expr = 1,
  track.progress = FALSE,
  EventType
)
}
\arguments{
\item{MarvelObject}{Marvel object. S3 object generated from \code{CreateMarvelObject} function.}

\item{min.cells}{Numeric value. The minimum number of cells in which the gene is expressed for the gene to included for splicing event detected and quantification. To be used in conjunction with \code{min.expr} argument. Default value is \code{50}.}

\item{min.expr}{Numeric value. The minimum expression value for the gene to be considered to be expressed in a cell. Default value is \code{1}.}

\item{track.progress}{Logical. If set to \code{TRUE}, progress bar will appear to track the progress of the rate-limiting step of this function, which is the extraction of the final exon-exon junctions. Default value is \code{FALSE}. Only applicable when \code{EventType} set to \code{"ALE"} or \code{"AFE"}.}

\item{EventType}{Character string. Indicate which splicing event type to calculate the PSI values for. Can take value \code{"ALE"} or \code{"AFE"}.}
}
\value{
An object of class S3 with new slot \code{MarvelObject$SpliceFeature$ALE} or \code{MarvelObject$SpliceFeature$AFE}.
}
\description{
Detects splicing events, specifically alternative first and last exons (AFE, ALE) from GTF. This is a wrapper function for \code{DetectEvents.ALE} and \code{DetectEvents.AFE} functions.
}
\examples{
marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))

marvel.demo <- DetectEvents(MarvelObject=marvel.demo,
                            min.cells=5,
                            min.expr=1,
                            track.progress=FALSE,
                            EventType="AFE"
                            )
}
