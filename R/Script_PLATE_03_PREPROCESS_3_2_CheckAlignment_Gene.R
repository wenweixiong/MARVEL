#' @title Check Alignment for gene data
#'
#' @description \code{CheckAlignment} checks if the metadata aligns with the columns and rows of the matrix for gene data.
#'
#' @details This function checks if the \code{sample.id} column of \code{MarvelObject$GenePheno} aligns with the column names of matrix in \code{MarvelObject$Exp}. To do this, the function subset overlapping sample IDs present in both the phenoData and matrix. This function also checks if the \code{gene_id} column of the \code{MarvelObject$GeneFeature} aligns with the \code{gene_id} column of matrix \code{MarvelObject$Exp}. To do this, the function subset overlappign gene IDs present in both the featureData and matrix.
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
#' marvel <- CheckAlignment.Exp(MarvelObject=marvel)

CheckAlignment.Exp <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    
    # Example argument
    #MarvelObject <- marvel
    
    #########################################################################
    ############################### PHENODATA ###############################
    #########################################################################
    
    # Print progress
    print("Checking phenoData...")

    # Retrieve phenoData
    df.pheno <- MarvelObject$GenePheno
    
    print(paste(length(df.pheno$sample.id), " samples identified in phenoData", sep=""))
    
    # Retrieve matrix
    df <- MarvelObject$Exp
    
    # Retrieve overlapping samples IDs
    sample.ids.phenoData <- df.pheno$sample.id
    sample.ids.matrix <- names(df)[-1]
    overlap <- intersect(sample.ids.phenoData, sample.ids.matrix)
        
    # Report progress
    print(paste(length(sample.ids.phenoData), " samples identified in phenoData ", sep=""))
    print(paste(length(sample.ids.matrix), " samples identified in matrix(s) ", sep=""))
        
    # Subset overlaps
    df.pheno <- df.pheno[which(df.pheno$sample.id %in% overlap), ]
    df <- df[,c("gene_id", overlap)]
    
    # Report progress
    print(paste(length(overlap), " overlapping samples retained ", sep=""))
        
    # Check alignment
    print("Checking alignment...")
    
    sample.ids.phenoData <- df.pheno$sample.id
    sample.ids.matrix <- names(df)[-1]
         
    index.l <- table(sample.ids.phenoData==sample.ids.matrix)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
    
        print("phenoData sample IDs and matrix column names MATCHED")
        
    } else {
        
        
        print("phenoData sample IDs and matrix column names NOT MATCHED")

    }
    
    # Save into MARVEL object
    MarvelObject$Exp <- df
    MarvelObject$GenePheno <- df.pheno

    #########################################################################
    ############################# FEATUREDATA ###############################
    #########################################################################
    
    # Print progress
    print("Checking featureData...")

    # Retrieve phenoData
    df.feature <- MarvelObject$GeneFeature
    
    print(paste(length(df.feature$gene_id), " genes identified in featureData", sep=""))
    
    # Retrieve matrix
    df <- MarvelObject$Exp
    
    # Retrieve overlapping genes
    gene.ids.featureData <- df.feature$gene_id
    gene.ids.matrix <- df$gene_id
    overlap <- intersect(gene.ids.featureData, gene.ids.matrix)
        
    # Report progress
    print(paste(length(gene.ids.featureData), " genes identified in featureData ", sep=""))
    print(paste(length(gene.ids.matrix), " genes identified in matrix(s) ", sep=""))
        
    # Subset overlaps
    df.feature <- df.feature[which(df.feature$gene_id %in% overlap), ]
    df <- df[which(df$gene_id %in% overlap),]
    
    # Report progress
    print(paste(length(overlap), " overlapping genes retained ", sep=""))
    
    # Check alignment
    print("Checking alignment...")
    
    gene.ids.featureData <- df.feature$gene_id
    gene.ids.matrix <- df$gene_id
         
    index.l <- table(gene.ids.featureData==gene.ids.matrix)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
    
        print("phenoData sample IDs and matrix column names MATCHED")
        
    } else {
        
        
        print("phenoData sample IDs and matrix column names NOT MATCHED")

    }

    # Save into MARVEL object
    MarvelObject$Exp <- df
    MarvelObject$GeneFeature <- df.feature
    
    #########################################################################
    
    # Return final object
    return(MarvelObject)
    
}


