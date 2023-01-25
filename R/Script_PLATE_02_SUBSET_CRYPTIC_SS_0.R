#' @title Subset cryptic alternative 5' and 3' splice site (A5/3SS)
#'
#' @description Subset A5SS and A3SS splicing events whose distance to the canonical are within a user-defined distance.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param DistanceToCanonical Numeric value. Distance between alternative and canonical splice site, below which, the alternative splice site will be defined as cryptic and retained.
#' @param EventType Character string. Type of splicing event. Only applicable to \code{"A5SS"} or \code{"A3SS"}.
#'
#' @return An object of class S3 with updated slot \code{$SpliceFeature$A5SS} or \code{$SpliceFeature$A3SS} and new slot \code{$DistanceToCanonical$A5SS} or \code{$DistanceToCanonical$A3SS}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export

SubsetCrypticSS <- function(MarvelObject, DistanceToCanonical=100, EventType) {


    if(EventType=="A5SS") {
        
        SubsetCrypticSS.A5SS(MarvelObject=MarvelObject,
                             DistanceToCanonical=DistanceToCanonical
                            )
        
    } else if(EventType=="A3SS") {
        
        SubsetCrypticSS.A3SS(MarvelObject=MarvelObject,
                             DistanceToCanonical=DistanceToCanonical
                            )
        
        
    }
   

}


