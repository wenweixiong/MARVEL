#' @title Compute percent spliced-in (PSI) values
#'
#' @description Validate splicing events and subsequently computes percent spliced-in (PSI) values these high-quality splicing events. This is a wrapper function for \code{ComputePSI.SE}, \code{ComputePSI.MXE}, \code{ComputePSI.A5SS}, \code{ComputePSI.A3SS}, \code{ComputePSI.RI}, \code{ComputePSI.AFE}, and \code{ComputePSI.ALE} functions.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param CoverageThreshold Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.
#' @param EventType Character string. Indicate which splicing event type to calculate the PSI values for. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param UnevenCoverageMultiplier Numeric value. Maximum allowable fold difference between two included junction counts for SE or two included or two excluded junction counts for MXE. Only applicable when \code{EventType} set to \code{"SE"} or \code{"MXE"}, respectively.
#' @param thread Numeric value. Only applicable when \code{EventType} set to \code{"RI"} Set number of threads..
#' @param read.length Numeric value. The length of read.Only applicable when \code{EventType} set to \code{"RI"}. This number will be specific to the sequencing mode. E.g. read length should be set to \code{150} when samples were sequenced in 150bp paired-end or single-end. This option should only be specified when users used read-counting approach for computing intron counts. The option should be left with its default value \code{1} when users tabulated the per-base count and summed them up to arrive at the intron counts.
#'
#' @return An object of class S3 with new slots \code{$SpliceFeatureValidated} and \code{$PSI}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- ComputePSI(MarvelObject=marvel.demo,
#'                           CoverageThreshold=10,
#'                           EventType="SE",
#'                           UnevenCoverageMultiplier=10
#'                           )

ComputePSI <- function(MarvelObject, CoverageThreshold, EventType, thread=NULL, UnevenCoverageMultiplier=10, read.length=1) {

    # Define arguments
    MarvelObject <- MarvelObject
    CoverageThreshold <- CoverageThreshold
    thread <- thread
    EventType <- EventType
    UnevenCoverageMultiplier <- UnevenCoverageMultiplier
    read.length <- read.length
    
    if(EventType=="SE") {
        
        ComputePSI.SE(MarvelObject=MarvelObject,
                      CoverageThreshold=CoverageThreshold,
                      UnevenCoverageMultiplier=UnevenCoverageMultiplier
                      )

    } else if(EventType=="MXE") {
        
        ComputePSI.MXE(MarvelObject=MarvelObject,
                       CoverageThreshold=CoverageThreshold,
                       UnevenCoverageMultiplier=UnevenCoverageMultiplier
                      )
                      
    } else if(EventType=="A5SS") {
        
        ComputePSI.A5SS(MarvelObject=MarvelObject,
                        CoverageThreshold=CoverageThreshold
                      )
    
    } else if(EventType=="A3SS") {
        
        ComputePSI.A3SS(MarvelObject=MarvelObject,
                        CoverageThreshold=CoverageThreshold
                      )
                      
    } else if(EventType=="RI") {
        
        ComputePSI.RI(MarvelObject=MarvelObject,
                        CoverageThreshold=CoverageThreshold,
                        IntronCounts=MarvelObject$IntronCounts,
                        thread=thread,
                        read.length=read.length
                        )
                      
    } else if(EventType=="ALE") {
        
        ComputePSI.ALE(MarvelObject=MarvelObject,
                       CoverageThreshold=CoverageThreshold
                       )
                       
    } else if(EventType=="AFE") {
        
        ComputePSI.AFE(MarvelObject=MarvelObject,
                       CoverageThreshold=CoverageThreshold
                       )
                       
    }
        
    
}
