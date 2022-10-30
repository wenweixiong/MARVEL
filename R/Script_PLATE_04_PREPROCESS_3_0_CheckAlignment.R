#' @title Pre-flight check
#'
#' @description Checks if the metadata aligns with the columns and rows of the matrix for splicing or gene data. This is a wrapper function for \code{CheckAlignment.PSI}, \code{CheckAlignment.Exp}, \code{CheckAlignment.PSI.Exp}, and \code{CheckAlignment.SJ}.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param level Character string. Indicate \code{"SJ"}, \code{"splicing"} or \code{"gene"} for splice junction, splicing or gene data, respectively. \code{"SJ"} typically specified before computing PSI values. \code{"splicing"} or \code{"gene"} typically specified after computing PSI values.
#'
#' @return An object of class S3 with updated slots \code{MarvelObject$SpliceJunction}, \code{MarvelObject$IntronCoverage}, \code{MarvelObject$SplicePheno}, \code{MarvelObject$SpliceFeatureValidated}, and \code{MarvelObject$PSI} or \code{MarvelObject$GenePheno}, \code{MarvelObject$GeneFeature}, and \code{MarvelObject$Gene} are updated for splicing or gene data, respectively.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- CheckAlignment(MarvelObject=marvel.demo,
#'                               level="SJ"
#'                               )

CheckAlignment <- function(MarvelObject, level) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    level <- level
        
    if(level=="SJ") {
        
        CheckAlignment.SJ(MarvelObject)
        
    } else if(level=="splicing") {
        
        CheckAlignment.PSI(MarvelObject)
        
    } else if(level=="gene") {
        
        CheckAlignment.Exp(MarvelObject)
        
    } else if(level=="splicing and gene")
    
        CheckAlignment.PSI.Exp(MarvelObject)
            
}


