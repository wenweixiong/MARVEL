#' @title Transform gene expression Values
#'
#' @description Transforms gene expression values and censor lowly-expressing genes.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment} function.
#' @param offset Numeric value. To indicate the value to add to the expression values before log transformation. The only option for this argument is 1.
#' @param transformation Character string. To indicate the type of transformation to use on the expression values after offsetting the values. The only option for this argument is \code{log2}.
#' @param threshold.lower Numeric value. To indicate the value below which the expression values will be censored, i.e. re-coded as 0, after offsetting and transforming the values. The only option for this argument is \code{1}.
#'
#' @return An object of class S3 with updated slot \code{MarvelObject$Exp}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- TransformExpValues(MarvelObject=marvel.demo,
#'                                   offset=1,
#'                                   transformation="log2",
#'                                   threshold.lower=1
#'                                   )

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
    
    # Check if data has been log2-transformed
    max.value <- max(unlist(df[,-1]))

    if(max.value < 20){
        
        message(paste("Maximum gene expression value detected as ", round(max.value, digits=2), ". It would seem that your gene expression values have been transformed. Hence, no transformation performed here.", sep="" ))
        
        return(MarvelObject)
        
    }
        
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
    message("Gene expression values offset by 1 and then log2-transformed. Transformed values below 1 re-coded as 0")
    
    # Save into MARVEL object
    MarvelObject$Exp <- df
    
    #########################################################################
    
    # Return final object
    return(MarvelObject)
    
}


