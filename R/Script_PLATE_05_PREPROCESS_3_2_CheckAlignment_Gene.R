#' @title Check gene data
#'
#' @description Checks if the metadata aligns with the columns and rows of the matrix for gene data.
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
#' marvel.demo <- CheckAlignment.Exp(MarvelObject=marvel.demo)

CheckAlignment.Exp <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- MarvelObject$GeneFeature
    df <- MarvelObject$Exp
    
    # Example argument
    #MarvelObject <- marvel
    #df.pheno <- MarvelObject$SplicePheno
    #df.feature <- MarvelObject$GeneFeature
    #df <- MarvelObject$Exp
    
    #########################################################################
    ############################### PHENODATA ###############################
    #########################################################################
        
    # Retrieve overlapping samples IDs
    sample.ids.phenoData <- df.pheno$sample.id
    sample.ids.matrix <- names(df)[-1]
    overlap <- intersect(sample.ids.phenoData, sample.ids.matrix)
        
    # Report progress
    message(paste(length(sample.ids.phenoData), " samples (cells) identified in sample metadata ", sep=""))
    message(paste(length(sample.ids.matrix), " samples (cells) identified in matrix ", sep=""))
        
    # Subset overlaps
    df.pheno <- df.pheno[which(df.pheno$sample.id %in% overlap), ]
    df <- df[,c("gene_id", overlap)]
    
    # Report progress
    message(paste(length(overlap), " overlapping samples (cells) retained ", sep=""))
        
    # Check alignment
    message("Checking alignment...")
    
    sample.ids.phenoData <- df.pheno$sample.id
    sample.ids.matrix <- names(df)[-1]
         
    index.l <- table(sample.ids.phenoData==sample.ids.matrix)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
    
        message("sample IDs in sample metadata and matrix column names MATCHED")
        
    } else {
        
        
        message("sample IDs in sample metadata and matrix column names NOT MATCHED")

    }

    #########################################################################
    ############################# FEATUREDATA ###############################
    #########################################################################
    
    # Print progress
    message(paste(length(df.feature$gene_id), " genes identified in gene metadata", sep=""))
    
    # Retrieve overlapping genes
    gene.ids.featureData <- df.feature$gene_id
    gene.ids.matrix <- df$gene_id
    overlap <- intersect(gene.ids.featureData, gene.ids.matrix)
        
    # Report progress
    message(paste(length(gene.ids.featureData), " genes identified in gene metadata", sep=""))
    message(paste(length(gene.ids.matrix), " genes identified in matrix", sep=""))
        
    # Subset overlaps
    df.feature <- df.feature[which(df.feature$gene_id %in% overlap), ]
    df <- df[which(df$gene_id %in% overlap),]
    
    # Report progress
    message(paste(length(overlap), " overlapping genes retained", sep=""))
    
    # Check alignment
    message("Checking alignment...")
    
    gene.ids.featureData <- df.feature$gene_id
    gene.ids.matrix <- df$gene_id
         
    index.l <- table(gene.ids.featureData==gene.ids.matrix)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
    
        message("gene IDs in gene metadata and matrix row names MATCHED")
        
    } else {
        
        
        message("gene IDs in gene metadata and matrix row names NOT MATCHED")

    }

    #########################################################################
    
    # Update slots
    MarvelObject$Exp <- df
    MarvelObject$SplicePheno <- df.pheno
    MarvelObject$GeneFeature <- df.feature
    
    # Return final object
    return(MarvelObject)
    
}


