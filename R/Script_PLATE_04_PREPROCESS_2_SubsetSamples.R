#' @title Subset samples (cells)
#'
#' @description Subsets specific samples (cells) from sample metadata.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param sample.ids Vector of character strings. Sample IDs to subset.
#'
#' @return An object of class S3 with updated slot \code{MarvelObject$SplicePheno}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' sample.ids <- sample(marvel.demo$SplicePheno$sample.id, size=10)
#'
#' marvel.demo <- SubsetSamples(MarvelObject=marvel.demo,
#'                              sample.ids=sample.ids
#'                              )

SubsetSamples <- function(MarvelObject, sample.ids) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df.pheno <- MarvelObject$SplicePheno
    sample.ids <- sample.ids
    
    # Example argument
    #MarvelObject <- marvel
    #df.pheno <- MarvelObject$SplicePheno
    #sample.id <- sample.ids
    
    # Track progress
    message(paste(length(df.pheno$sample.id), " samples identified in sample metadata", sep=""))
    
    # Subset overlapping samples in sample metadata
    overlap <- intersect(df.pheno$sample.id, sample.ids)
    df.pheno <- df.pheno[which(df.pheno$sample.id %in% overlap), ]
    
    # Track progress
    message(paste(length(df.pheno$sample.id), " samples retained in sample metadata", sep=""))
    
    #########################################################################
    
    # Update slots
    MarvelObject$SplicePheno <- df.pheno
    
    # Return final object
    return(MarvelObject)
    
}


