#' @title Check splice junction data
#'
#' @description Checks if the metadata aligns with the columns and rows of the matrix for splice junction data prior to PSI computation.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#'
#' @return An object of class S3 with updated slots \code{MarvelObject$SplicePheno}, \code{MarvelObject$PSI} and \code{MarvelObject$IntronCounts}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- CheckAlignment.SJ(MarvelObject=marvel.demo)

CheckAlignment.SJ <- function(MarvelObject) {

    # Define arguments
    MarvelObject <- MarvelObject
        
    # Check alignment: SJ matrix and phenoData
        # Retrieve data
        df.pheno <- MarvelObject$SplicePheno
        df <- MarvelObject$SpliceJunction
        
        # Report overlapping samples
        overlap <- intersect(df.pheno$sample.id, names(df)[-1])
        message(paste(length(df.pheno$sample.id), " samples (cells) identified in SJ phenoData", sep=""))
        message(paste(length(names(df)[-1]), " samples (cells) identified in SJ count matrix", sep=""))
        message(paste(length(overlap), " overlapping samples (cells) identified", sep=""))
        
        # Subset overlapping samples
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% overlap), ]
        df <- df[, c("coord.intron", df.pheno$sample.id)]
        
        # Check alignment
        index.l <- table(df.pheno$sample.id==names(df)[-1])
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
        
        if(index.true==1 & index.false==0) {
            
           message("sample IDs in sample metadata and SJ count matrix column names MATCHED")
           
        } else {
           
           
           message("sample IDs in sample metadata and SJ count matrix column names NOT MATCHED")
           
        }
        
        
    # Check alignment: Additional check for RI
    if(inherits(MarvelObject$IntronCounts, "data.frame", TRUE)==1) {
        
        # Report progress
        message("Additional checks for intron count matrix...")
        
        # Retrieve data
        df.intron.counts <- MarvelObject$IntronCounts
        
        # Report overlapping samples
        overlap <- intersect(df.pheno$sample.id, names(df.intron.counts)[-1])
        message(paste(length(df.pheno$sample.id), " samples (cells) identified in SJ phenoData", sep="" ))
        message(paste(length(names(df.intron.counts)[-1]), " samples (cells) identified in intron count matrix", sep=""))
        message(paste(length(overlap), " overlapping samples (cells) identified", sep=""))
        
        # Subset overlapping samples
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% overlap), ]
        df.intron.counts <- df.intron.counts[, c("coord.intron", df.pheno$sample.id)]
        df <- df[, c("coord.intron", df.pheno$sample.id)]
        
        # Check alignment
        index.l <- table(df.pheno$sample.id==names(df.intron.counts)[-1])
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
        
        if(index.true==1 & index.false==0) {
            
           message("sample IDs in sample metadata and SJ count matrix and intron count column names MATCHED")
           
        } else {
           
           
           message("sample IDs in sample metadata and SJ count matrix and intron count matrix column names NOT MATCHED")
           
        }
                
    }
    
    ##################################################
    
    # Update slots
    MarvelObject$SplicePheno <- df.pheno
    MarvelObject$SpliceJunction <- df
    MarvelObject$IntronCounts <- df.intron.counts
    
    # Return final object
    return(MarvelObject)
    
}
