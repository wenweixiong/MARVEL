#' @title Check splicing and gene data against each other
#'
#' @description Subsets overlapping samples between splicing and gene data.
#'
#' @param MarvelObject S3 object generated from \code{CheckAlignment.PSI} and \code{CheckAlignment.Exp} function.
#'
#' @return An object of class S3 with updated slots \code{MarvelObject$SplicePheno}, \code{MarvelObject$PSI}, \code{MarvelObject$GenePheno}, and \code{MarvelObject$Exp}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- CheckAlignment.PSI.Exp(MarvelObject=marvel.demo)

CheckAlignment.PSI.Exp <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df.pheno <- MarvelObject$SplicePheno
    df.exp <- MarvelObject$Exp
    
    # Example argument
    #MarvelObject <- marvel
    #df.pheno <- MarvelObject$SplicePheno
    #df.exp <- MarvelObject$Exp
    
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
    
    # Retrieve overlapping samples IDs
        # Gene expression matrix
        sample.ids.exp <- names(df.exp)[-1]
    
        # Matrix
        .list <- list()
        
        for(i in 1:length(event.types)) {
             
            .list[[i]] <- names(MarvelObject$PSI[[event.types[i]]])[-1]
            
        }
        
        sample.ids.psi <- unique(unlist(.list))
        
        # Overlaps
        overlap <- intersect(sample.ids.exp, sample.ids.psi)
        
        # Report progress
        message(paste(length(sample.ids.psi), " samples (cells) identified in psi matrix(s)", sep=""))
        message(paste(length(sample.ids.exp), " samples (cells) identified in gene expression matrix", sep=""))
        message(paste(length(overlap), " overlapping samples (cells) identified and retained", sep=""))
        
    # Subset overlapping samples IDs
        # PSI matrix
        for(i in 1:length(event.types)) {
             
             df <- MarvelObject$PSI[[event.types[i]]]
            
             df.small <- df[, c("tran_id", overlap)]
             
             MarvelObject$PSI[[event.types[i]]] <- df.small
            
        }
        
        # Gene matrix
        df.exp <- df.exp[,c("gene_id", overlap)]
    
        
    # Check alignment
    for(i in 1:length(MarvelObject$PSI)) {
        
        # Retrieve PSI matrix
        df.psi <- MarvelObject$PSI[[i]]
        
        index.l <- table(names(df.psi)[-1]==names(df.exp)[-1])
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message(paste("sample IDs in columns of gene matrix and ", names(MarvelObject$PSI)[i], " PSI matrix MATCHED", sep=""))
            
        } else {
            
            message(paste("sample IDs in columns of gene matrix and ", names(MarvelObject$PSI)[i], " PSI matrix NOT MATCHED", sep=""))

        }
        
    }
    
    #########################################################################
    
    # Update slots
    MarvelObject$SplicePheno <- df.pheno
    MarvelObject$MarvelObject$Exp <- df.exp
    
    # Return final object
    return(MarvelObject)
    
}


