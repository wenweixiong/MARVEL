#' @title Subset Samples (Cells)
#'
#' @description \code{SubsetSamples} subsets samples (cells) based on user-defined criteria for splicing or gene data.
#'
#' @details This function subset samples (cells) based on user-defined criteria, e.g. samples passing sequencing QC. These criteria are reflected in the columns of the \code{MarvelObject$SplicePheno} or \code{MarvelObject$GenePheno} slot.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} function.
#' @param columns Character string. To indicate which columns in the \code{MarvelObject$SplicePheno} or \code{MarvelObject$GenePheno} slot to refer to when filtering samples (cells).
#' @param variables List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{columns} argument.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for splicing or gene data, respectively.
#'
#' @export
#' @return An object of class S3. The original \code{MarvelObject$SplicePheno} or \code{MarvelObject$GenePheno} slot is replaced with the new filtered phenoData for splicing or gene data, respectively.
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
#' marvel <- SubsetSamples(MarvelObject=marvel,
#'                         columns=c("qc.seq", "sample.type", "cell.type"),
#'                         variables=list("pass", "Single Cell", c("iPSC", "Endoderm")),
#'                         level="splicing"
#'                         )

SubsetSamples <- function(MarvelObject, columns, variables, level) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    columns <- columns
    variables <- variables
    level <- level
        
    if(level=="splicing") {
        
        SubsetSamples.PSI(MarvelObject, columns, variables)
        
    } else if(level=="gene") {
        
        SubsetSamples.Gene(MarvelObject, columns, variables)
        
    }

            
}


