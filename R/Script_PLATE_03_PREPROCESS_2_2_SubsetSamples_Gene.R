#' @title Subset Samples (Cells) for gene data
#'
#' @description \code{SubsetSamples} subsets samples (cells) based on user-defined criteria for gene data.
#'
#' @details This function subset samples (cells) based on user-defined criteria, e.g. samples passing sequencing QC. These criteria are reflected in the columns of the \code{MarvelObject$GenePheno} slot.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} function.
#' @param columns Character string. To indicate which columns in the \code{MarvelObject$GenePheno} slot to refer to when filtering samples (cells).
#' @param variables List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{columns} argument.
#'
#' @export
#'
#' @return An object of class S3. The original \code{MarvelObject$GenePheno} slot is replaced with the new filtered phenoData.
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
#' marvel <- SubsetSamples.Gene(MarvelObject=marvel,
#'                              columns=c("qc.seq", "sample.type", "cell.type"),
#'                              variables=list("pass", "Single Cell", c("iPSC", "Endoderm"))
#'                              )

SubsetSamples.Gene <- function(MarvelObject, columns, variables) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    columns <- columns
    variables <- variables
    
    # Example argument
    #MarvelObject <- marvel
    #columns <- c("qc.seq", "control.well", "sample.type")
    #variables <- c("pass", "FALSE", "Single Cell")
        
    # Retrieve phenoData
    df.pheno <- MarvelObject$GenePheno
    
    # Track progress
    print(paste(length(df.pheno$sample.id), " samples identified in phenoData", sep=""))
    
    # Retrieve relevant cell types
    index.list <- list()
    
    for(i in 1:length(columns)) {
        
        index.list[[i]] <- which(df.pheno[[columns[i]]] %in% variables[[i]])

    }
   
    index <- Reduce(intersect, index.list)
    
    # Subset relevant cell types and save into MARVEL object: phenoData
    df.pheno <- df.pheno[index, ]
    MarvelObject$GenePheno <- df.pheno
    
    # Track progress
    print(paste(length(df.pheno$sample.id), " samples retained in phenoData", sep=""))
    
    #########################################################################
    
    # Return final object
    return(MarvelObject)
    
}


