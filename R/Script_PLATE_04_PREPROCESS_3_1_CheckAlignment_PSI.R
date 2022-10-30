#' @title Check splicing data
#'
#' @description Checks if the metadata aligns with the columns and rows of the matrix for splicing data.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#'
#' @return An object of class S3 with updated slots \code{MarvelObject$SplicePheno}, \code{MarvelObject$SpliceFeature}, and \code{MarvelObject$PSI}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- CheckAlignment.PSI(MarvelObject=marvel.demo)

CheckAlignment.PSI <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df.pheno <- MarvelObject$SplicePheno
    
    # Example argument
    #MarvelObject <- marvel
    #df.pheno <- MarvelObject$SplicePheno
    
    # Retrieve non-empty matrix event types
    ncol <- c(ncol(MarvelObject$PSI[["SE"]]),
              ncol(MarvelObject$PSI[["MXE"]]),
              ncol(MarvelObject$PSI[["RI"]]),
              ncol(MarvelObject$PSI[["A5SS"]]),
              ncol(MarvelObject$PSI[["A3SS"]]),
              ncol(MarvelObject$PSI[["ALE"]]),
              ncol(MarvelObject$PSI[["AFE"]])
              )
    event.types <- c("SE", "MXE", "RI", "A5SS", "A3SS", "ALE", "AFE")
    event.types <- event.types[which(ncol >= 2)]
    
    #########################################################################
    ############################### PHENODATA ###############################
    #########################################################################
    
    # Track progress
    message(paste(length(df.pheno$sample.id), " samples (cells) identified in sample metadata", sep=""))
    
    # Retrieve overlapping samples IDs
        # phenoData
        sample.ids.phenoData <- df.pheno$sample.id
    
        # Matrix
        .list <- list()
        
        for(i in 1:length(event.types)) {
             
            .list[[i]] <- names(MarvelObject$PSI[[event.types[i]]])
            
            
        }
        
        sample.ids.matrix <- unique(unlist(.list))
        sample.ids.matrix <- sample.ids.matrix[which(sample.ids.matrix != "tran_id")]
        
        # Overlaps
        overlap <- intersect(sample.ids.phenoData, sample.ids.matrix)
        
        # Report progress
        message(paste(length(sample.ids.phenoData), " samples (cells) identified in sample metadata", sep=""))
        message(paste(length(sample.ids.matrix), " samples (cells) identified in matrix(s) ", sep=""))
        
    # Subset and return MARVEL object
        # phenoData
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% overlap), ]
                
        # Matrix
        for(i in 1:length(event.types)) {
             
             df <- MarvelObject$PSI[[event.types[i]]]
            
             df.small <- df[, c("tran_id", overlap)]
             
             MarvelObject$PSI[[event.types[i]]] <- df.small
            
        }
        
        # Report progress
        message(paste(length(overlap), " overlapping samples (cells) retained", sep=""))
        
    # Check alignment
    message("Checking alignment...")
    
    sample.ids.phenoData <- df.pheno$sample.id
    
    for(i in 1:length(event.types)) {
         
         sample.ids.matrix <- names(MarvelObject$PSI[[event.types[i]]])[-1]
         
         index.l <- table(sample.ids.phenoData==sample.ids.matrix)
         index.true <- length(which(names(index.l)==TRUE))
         index.false <- length(which(names(index.l)==FALSE))
         
         if(index.true==1 & index.false==0) {
             
            message(paste("sample IDs in sample metadata and matrix column names MATCHED for ", event.types[i], sep=""))
            
         } else {
            
            
            message(paste("sample IDs in sample metadata and matrix column names NOT MATCHED for ", event.types[i], sep=""))
            
         }
        
    }
            
    #########################################################################
    ############################# FEATUREDATA ###############################
    #########################################################################
    
    # Subset and return MARVEL object
    for(i in 1:length(event.types)) {
        
        # Report progress
        message(paste("Checking for ", event.types[i], "...", sep=""))
        
        # Retrieve matrix
        df <- MarvelObject$PSI[[event.types[i]]]
        
        # Retrieve featureData
        df.feature <- MarvelObject$SpliceFeatureValidated[[event.types[i]]]
        
        # Retrieve overlaps
        tran.ids.featureData <- df.feature$tran_id
        tran.ids.matrix <- df$tran_id
        overlap <- intersect(tran.ids.featureData, tran.ids.matrix )
        
        # Report progress
        message(paste(length(tran.ids.featureData), " transcripts identified in transcript metadata", sep=""))
        message(paste(length(tran.ids.matrix), " transcripts identified in matrix", sep=""))
        
        # Subset overlaps
        df.feature <- df.feature[which(df.feature$tran_id %in% overlap), , drop=FALSE]
        df <- df[which(df$tran_id %in% overlap), ]
        
        # Report progress
        message(paste(length(overlap), " overlapping transcripts retained", sep=""))
        
        # Save into MARVEL object
        MarvelObject$SpliceFeatureValidated[[event.types[i]]] <- df.feature
        MarvelObject$PSI[[event.types[i]]] <- df
                    
    }
    
    # Check alignment
    message("Checking alignment...")
    
    for(i in 1:length(event.types)) {
        
         # featureData
         tran.ids.featureData <- MarvelObject$SpliceFeatureValidated[[event.types[i]]]$tran_id

         # Matrix
         tran.ids.matrix <- MarvelObject$PSI[[event.types[i]]]$tran_id
         
         index.l <- table(tran.ids.featureData==tran.ids.matrix)
         index.true <- length(which(names(index.l)==TRUE))
         index.false <- length(which(names(index.l)==FALSE))
         
         if(index.true==1 & index.false==0) {
             
            message(paste("Transcript IDs in transcript metadata and matrix row names MATCHED for ", event.types[i], sep=""))
            
         } else {
            
            
            message(paste("Transcript IDs in transcript metadata and matrix row names NOT MATCHED for ", event.types[i], sep=""))
            
         }
        
    }

    #########################################################################
    
    # Return final object
    return(MarvelObject)
    
}


