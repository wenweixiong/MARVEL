#' @title Compute Percent Spliced-in (PSI) Values
#'
#' @description \code{ComputePSI} computes percent spliced-in (PSI) values for a specified splicing event type.
#'
#' @details This function computes the PSI values for a splicing event as specified in the \code{EventType} argument. Splicing events provided in \code{SpliceFeature} data frame will first be cross-checked against the splice junctions provided in \code{SpliceJunction} data frame. Only events whose junctions are found in \code{SpliceJunction} are retained. The formula for computing PSI is the number of junction reads supporting the included isoform divided by the total number of reads supporting both included and excluded isoforms.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} function.
#' @param CoverageThreshold Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.
#' @param EventType Character string. Indicate which splicing event type to calculate the PSI values for. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param IntronCountsFile Data frame containing per base coverage of introns. Only applicable when \code{EventType} set to \code{"RI"}. First column should be named \code{coord.intron} and indicate the per base intron position in the form of chr:position. Subsequent columns should contain the per base coverage for each sample. These counts can be deteted using external softwares such as Bedtools etc..
#' @param thread Numeric value. Set number of threads.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to two new slots. \code{$SpliceFeatureValidated} contains the validated splicing event metadata. \code{$PSI} contains the computed PSI values for the validated splicing events.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @import methods
#' @import parallel
#' @importFrom plyr join
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' marvel <- ComputePSI(MarvelObject=marvel,
#'                      CoverageThreshold=10,
#'                     EventType="SE"
#'                     )
#'
#' marvel$SpliceFeatureValidated$SE
#' marvel$PSI$SE[,1:5]

ComputePSI <- function(MarvelObject, CoverageThreshold, EventType, IntronCountsFile=NULL, thread=NULL) {

    # Define arguments
    MarvelObject <- MarvelObject
    CoverageThreshold <- CoverageThreshold
    IntronCountsFile <- IntronCountsFile
    thread <- thread
    EventType <- EventType
    
    if(EventType=="SE") {
        
        ComputePSI.SE(MarvelObject=MarvelObject,
                      CoverageThreshold=CoverageThreshold
                      )

    } else if(EventType=="MXE") {
        
        ComputePSI.MXE(MarvelObject=MarvelObject,
                       CoverageThreshold=CoverageThreshold
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
                        IntronCountsFile=IntronCountsFile,
                        thread=thread
                        )
                      
    }
    
}
