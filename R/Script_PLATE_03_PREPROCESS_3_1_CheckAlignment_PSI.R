#' @title Check Alignment for splicing data
#'
#' @description \code{CheckAlignment} checks if the metadata aligns with the columns and rows of the matrix for splicing data.
#'
#' @details This function checks if the \code{sample.id} column of \code{MarvelObject$SplicePheno} aligns with the column names of matrix in \code{MarvelObject$PSI}. To do this, the function subset overlapping sample IDs present in both the phenoData and matrix. This function also checks if the \code{tran_id} column of the \code{MarvelObject$SpliceFeature} aligns with the \code{tran_id} column of matrix in \code{MarvelObject$PSI}. To do this, the function subset overlappign gene IDs present in both the featureData and matrix.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{SubsetSamples} function.
#'
#' @export
#'
#' @return An object of class S3. The original \code{MarvelObject$SplicePheno}, \code{MarvelObject$SpliceFeature}, and \code{MarvelObject$PSI} are updated
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
#' marvel <- CheckAlignment.PSI(MarvelObject=marvel)


CheckAlignment.PSI <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    
    # Example argument
    #MarvelObject <- marvel
    
    # Retrieve non-empty matrix event types
    ncol <- c(ncol(MarvelObject$PSI[["SE"]]),
              ncol(MarvelObject$PSI[["MXE"]]),
              ncol(MarvelObject$PSI[["RI"]]),
              ncol(MarvelObject$PSI[["A5SS"]]),
              ncol(MarvelObject$PSI[["A3SS"]])
              )
    event.types <- c("SE", "MXE", "RI", "A5SS", "A3SS")
    event.types <- event.types[which(ncol >= 2)]
    
    #########################################################################
    ############################### PHENODATA ###############################
    #########################################################################
    
    # Print progress
    print("Checking phenoData...")

    # Retrieve phenoData
    df.pheno <- MarvelObject$SplicePheno
    
    print(paste(length(df.pheno$sample.id), " samples identified in phenoData", sep=""))
    
    # Retrieve overlapping samples IDs
        # phenoData
        sample.ids.phenoData <- df.pheno$sample.id
    
        # Matrix
        .list <- list()
        
        for(i in 1:length(event.types)) {
             
            .list[[i]] <- names(MarvelObject$PSI[[event.types[i]]])
            
            
        }
        
        sample.ids.matrix <- unique(unlist(.list))
        
        # Overlaps
        overlap <- intersect(sample.ids.phenoData, sample.ids.matrix)
        
        # Report progress
        print(paste(length(sample.ids.phenoData), " samples identified in phenoData ", sep=""))
        print(paste(length(sample.ids.matrix), " samples identified in matrix(s) ", sep=""))
        
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
        print(paste(length(overlap), " overlapping samples retained ", sep=""))
        
    # Check alignment
    print("Checking alignment...")
    
    sample.ids.phenoData <- df.pheno$sample.id
    
    for(i in 1:length(event.types)) {
         
         sample.ids.matrix <- names(MarvelObject$PSI[[event.types[i]]])[-1]
         
         index.l <- table(sample.ids.phenoData==sample.ids.matrix)
         index.true <- length(which(names(index.l)==TRUE))
         index.false <- length(which(names(index.l)==FALSE))
         
         if(index.true==1 & index.false==0) {
             
            print(paste("phenoData sample IDs and matrix column names MATCHED for ", event.types[i], sep=""))
            
         } else {
            
            
            print(paste("phenoData sample IDs and matrix column names NOT MATCHED for ", event.types[i], sep=""))
            
         }
        
    }
            
    #########################################################################
    ############################# FEATUREDATA ###############################
    #########################################################################
    
    # Subset and return MARVEL object
    for(i in 1:length(event.types)) {
        
        # Report progress
        print(paste("Checking for ", event.types[i], sep=""))
        
        # Retrieve matrix
        df <- MarvelObject$PSI[[event.types[i]]]
        
        # Retrieve featureData
        df.feature <- MarvelObject$SpliceFeatureValidated[[event.types[i]]]
        
        # Retrieve overlaps
        tran.ids.featureData <- df.feature$tran_id
        tran.ids.matrix <- df$tran_id
        overlap <- intersect(tran.ids.featureData, tran.ids.matrix )
        
        # Report progress
        print(paste(length(tran.ids.featureData), " transcripts identified in featureData ", sep=""))
        print(paste(length(tran.ids.matrix), " transcripts identified in matrix", sep=""))
        
        # Subset overlaps
        df.feature <- df.feature[which(df.feature$tran_id %in% overlap), , drop=FALSE]
        df <- df[which(df$tran_id %in% overlap), ]
        
        # Report progress
        print(paste(length(overlap), " overlapping transcripts retained", sep=""))
        
        # Save into MARVEL object
        MarvelObject$SpliceFeatureValidated[[event.types[i]]] <- df.feature
        MarvelObject$PSI[[event.types[i]]] <- df
                    
    }
    
    # Check alignment
    print("Checking alignment...")
    
    for(i in 1:length(event.types)) {
        
         # featureData
         tran.ids.featureData <- MarvelObject$SpliceFeatureValidated[[event.types[i]]]$tran_id

         # Matrix
         tran.ids.matrix <- MarvelObject$PSI[[event.types[i]]]$tran_id
         
         index.l <- table(tran.ids.featureData==tran.ids.matrix)
         index.true <- length(which(names(index.l)==TRUE))
         index.false <- length(which(names(index.l)==FALSE))
         
         if(index.true==1 & index.false==0) {
             
            print(paste("featureData transcript IDs and matrix row names MATCHED for ", event.types[i], sep=""))
            
         } else {
            
            
            print(paste("featureData transcript IDs and matrix row names NOT MATCHED for ", event.types[i], sep=""))
            
         }
        
    }

    #########################################################################
    
    # Return final object
    return(MarvelObject)
    
}


