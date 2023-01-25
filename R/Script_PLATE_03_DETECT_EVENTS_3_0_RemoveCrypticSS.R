#' @title Remove cryptic alternative 5' and 3' splice site (A5/3SS)
#'
#' @description Subset A5SS and A3SS splicing events from alternative first and last exon (AFE, ALE), respectively, whose distance to the canonical are within a user-defined distance.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param DistanceToCanonical Numeric value. Distance between alternative and canonical splice site, below which, the alternative splice site will be defined as cryptic and removed.
#' @param EventType Character string. Type of splicing event. Only applicable to \code{"A5SS"} or \code{"A3SS"}.
#'
#' @return An object of class S3 with updated slot \code{$SpliceFeature$AFE} or \code{$SpliceFeature$ALE} and new slot \code{$DistanceToCanonical$AFE} or \code{$DistanceToCanonical$ALE}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export

RemoveCrypticSS <- function(MarvelObject, DistanceToCanonical=100, EventType) {

    if(EventType=="AFE") {
        
        RemoveCrypticSS.AFE(MarvelObject=MarvelObject,
                             DistanceToCanonical=DistanceToCanonical
                            )
        
    } else if(EventType=="ALE") {
        
        RemoveCrypticSS.ALE(MarvelObject=MarvelObject,
                             DistanceToCanonical=DistanceToCanonical
                            )
        
    }

}


