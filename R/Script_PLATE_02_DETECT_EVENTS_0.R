#' @title Detect Splicing Events
#'
#' @description Detects splicing events, specifically alternative first and last exons (AFE, ALE) from GTF. This is a wrapper function for \code{DetectEvents.ALE} and \code{DetectEvents.AFE} functions.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param min.cells Numeric value. The minimum number of cells in which the gene is expressed for the gene to included for splicing event detected and quantification. To be used in conjunction with \code{min.expr} argument. Default value is \code{50}.
#' @param min.expr Numeric value. The minimum expression value for the gene to be considered to be expressed in a cell. Default value is \code{1}.
#' @param track.progress Logical. If set to \code{TRUE}, progress bar will appear to track the progress of the rate-limiting step of this function, which is the extraction of the final exon-exon junctions. Default value is \code{FALSE}. Only applicable when \code{EventType} set to \code{"ALE"} or \code{"AFE"}.
#' @param EventType Character string. Indicate which splicing event type to calculate the PSI values for. Can take value \code{"ALE"} or \code{"AFE"}.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$SpliceFeature$ALE} or \code{MarvelObject$SpliceFeature$AFE}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- DetectEvents(MarvelObject=marvel.demo,
#'                             min.cells=5,
#'                             min.expr=1,
#'                             track.progress=FALSE,
#'                             EventType="AFE"
#'                             )

DetectEvents <- function(MarvelObject, min.cells=50, min.expr=1, track.progress=FALSE, EventType) {

    # Define arguments
    min.cells <- min.cells
    min.expr <- min.expr
    track.progress <- track.progress
    EventType <- EventType
    
    if(EventType=="ALE") {
        
        DetectEvents.ALE(MarvelObject=MarvelObject,
                         min.cells=min.cells,
                         min.expr=min.expr,
                         track.progress=track.progress
                         )

    } else if(EventType=="AFE") {
        
        DetectEvents.AFE(MarvelObject=MarvelObject,
                         min.cells=min.cells,
                         min.expr=min.expr,
                         track.progress=track.progress
                         )
                      
    }
        
    
}
