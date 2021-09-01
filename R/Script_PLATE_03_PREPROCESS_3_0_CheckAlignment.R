#' @title Check Alignment
#'
#' @description \code{CheckAlignment} checks if the metadata aligns with the columns and rows of the matrix for splicing or gene data.
#'
#' @details This function checks if the \code{sample.id} column of \code{MarvelObject$SplicePheno} or \code{MarvelObject$GenePheno} aligns with the column names of matrix in \code{MarvelObject$PSI} or \code{MarvelObject$Gene} for splicing or gene data, respectively. To do this, the function subset overlapping sample IDs present in both the phenoData and matrix. This function also checks if the \code{tran_id} or \code{gene_id} column of the \code{MarvelObject$SpliceFeature} or \code{MarvelObject$GeneFeature} aligns with the \code{tran_id} or \code{gene_id} column of matrix in \code{MarvelObject$PSI} or \code{MarvelObject$Gene} for splicing or gene data, respectively. To do this, the function subset overlappign gene IDs present in both the featureData and matrix.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{SubsetSamples} function.
#' @param level Character string. Indicate \code{"SJ"}, \code{"splicing"} or \code{"gene"} for splice junction, splicing or gene data, respectively. \code{"SJ"} typically specified before computing PSI values. \code{"splicing"} or \code{"gene"} typically specified after computing PSI values.
#' 
#' @export
#'
#' @return An object of class S3. The original \code{MarvelObject$SpliceJunction}, \code{MarvelObject$IntronCoverage}, \code{MarvelObject$SplicePheno}, \code{MarvelObject$SpliceFeatureValidated}, and \code{MarvelObject$PSI} or \code{MarvelObject$GenePheno}, \code{MarvelObject$GeneFeature}, and \code{MarvelObject$Gene} are updated for splicing or gene data, respectively.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @import methods
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- CheckAlignment(MarvelObject=marvel, level="SJ")


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


