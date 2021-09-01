#' @title Check Alignment for splice junction data
#'
#' @description \code{CheckAlignment} checks if the metadata aligns with the columns and rows of the matrix for splice junction data.
#'
#' @details This function checks if the \code{sample.id} column of \code{MarvelObject$SplicePheno} aligns with the column names of splice junction count matrix in \code{MarvelObject$SpliceJunction}. To do this, the function subset overlapping sample IDs present in both the phenoData and matrix. If the intron coverage matrix is found in \code{MarvelObject$IntronCounts}, additional checking will be done to check if the \code{sample.id} column of \code{MarvelObject$SplicePheno} and column ames of the splice junction count matrix in \code{MarvelObject$SpliceJunction} align with the column names of intron coverage matrix.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{SubsetSamples} function.
#'
#' @export
#'
#' @return An object of class S3. The original \code{MarvelObject$SplicePheno}, and \code{MarvelObject$PSI} are updated. \code{MarvelObject$IntronCounts} updated if intron coverage matrix detected.
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
#' marvel <- CheckAlignment.SJ(MarvelObject=marvel)

CheckAlignment.SJ <- function(MarvelObject) {

    # Define arguments
    MarvelObject <- MarvelObject
        
    # Check alignment: SJ matrix and phenoData
        # Retrieve data
        df.pheno <- MarvelObject$SplicePheno
        df <- MarvelObject$SpliceJunction
        
        # Report overlapping samples
        overlap <- intersect(df.pheno$sample.id, names(df)[-1])
        print(paste(length(df.pheno$sample.id), " samples (cells) identified in SJ phenoData", sep=""))
        print(paste(length(names(df)[-1]), " samples (cells) identified in SJ count matrix", sep=""))
        print(paste(length(overlap), " overlapping samples (cells) identified", sep=""))
        
        # Subset overlapping samples
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% overlap), ]
        df <- df[, c("coord.intron", df.pheno$sample.id)]
        
        # Check alignment
        index.l <- table(df.pheno$sample.id==names(df)[-1])
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
        
        if(index.true==1 & index.false==0) {
            
           print("phenoData SJ IDs and SJ count matrix column names MATCHED")
           
        } else {
           
           
           print("phenoData sample IDs and SJ count matrix column names NOT MATCHED")
           
        }
        
        
    # Check alignment: Additional check for RI
    if(class(MarvelObject$IntronCounts) == "data.frame") {
        
        # Report progress
        print("Additional checks for intron count matrix...")
        
        # Retrieve data
        df.intron.counts <- MarvelObject$IntronCounts
        
        # Report overlapping samples
        overlap <- intersect(df.pheno$sample.id, names(df.intron.counts)[-1])
        print(paste(length(df.pheno$sample.id), " samples (cells) identified in SJ phenoData", sep="" ))
        print(paste(length(names(df.intron.counts)[-1]), " samples (cells) identified in intron count matrix", sep=""))
        print(paste(length(overlap), " overlapping samples (cells) identified", sep=""))
        
        # Subset overlapping samples
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% overlap), ]
        df.intron.counts <- df.intron.counts[, c("coord.intron", df.pheno$sample.id)]
        df <- df[, c("coord.intron", df.pheno$sample.id)]
        
        # Check alignment
        index.l <- table(df.pheno$sample.id==names(df.intron.counts)[-1])
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
        
        if(index.true==1 & index.false==0) {
            
           print("phenoData SJ IDs and SJ count matrix and intron count column names MATCHED")
           
        } else {
           
           
           print("phenoData SJ IDs and SJ count matrix and intron count matrix column names NOT MATCHED")
           
        }
                
    }
    
    # Save into slots
    MarvelObject$SplicePheno <- df.pheno
    MarvelObject$SpliceJunction <- df
    MarvelObject$IntronCounts <- df.intron.counts
    
    # Return final object
    return(MarvelObject)
    
}
