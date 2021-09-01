#' @title Transform Gene Expression Values
#'
#' @description \code{TransformExpValues} transforms gene expression values.
#'
#' @details This function transforms gene expression values prior to downstream analysis, such as differential expression analysis. The values provided to this function are assumed to have already been normalised, e.g. normalised by per-cell library size and gene length and then multiplied by a scale factor in the case of transcript per million (TPM).
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject}, \code{SubsetSamples} or \code{CheckAlignment} function.
#' @param offset Numeric value. To indicate the value to add to the expression values before log transformation. The only option for this argument is 1.
#' @param transformation Character string. To indicate the type of transformation to use on the expression values after offsetting the values. The only option for this argument is \code{log2}.
#' @param threshold.lower Numeric value. To indicate the value below which the expression values will be censored, i.e. re-coded as 0, after offsetting and transforming the values. The only option for this argument is \code{1}.
#'
#' @export
#' @return An object of class S3. The original matirx in \code{MarvelObject$Exp} is replaced with the new offset-ed, transformed, and censored values.
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
#' marvel <- TransformExpValues(MarvelObject=marvel)

TransformExpValues <- function(MarvelObject, offset=1, transformation="log2", threshold.lower=1) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    offset <- offset
    transformation <- transformation
    threshold.lower <- threshold.lower
    
    # Example argument
    #MarvelObject <- marvel
    #offset <- 1
    #transformation <- "log2"
    #threshold.lower <- 1

    # Retrieve matrix
    df <- MarvelObject$Exp
    
    # Apply offset
    df[,-1] <- df[,-1] + offset
    
    # Apply log-transformation
    if(transformation=="log2") {
        
        df[,-1] <- log2(df[,-1])
                
    } else if(transformation=="log10") {
        
        df[,-1] <- log10(df[,-1])
        
    }
    
    # Recode values below lower threshold to 0
    df[,-1][df[,-1] < threshold.lower] <- 0
    
    # Track progress
    print("Gene expression values offset by 1 and then log2-transformed. Transformed values below 1 re-coded as 0")
    
    # Save into MARVEL object
    MarvelObject$Exp <- df
    
    #########################################################################
    
    # Return final object
    return(MarvelObject)
    
}


