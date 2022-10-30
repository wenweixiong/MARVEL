#' @title Differential gene expression analysis of specified gene
#'
#' @description Performs differential gene expression analysis specified gene across for all possible pairs of cell groups. The gene and cell groups were defined earlier in \code{adhocGene.TabulateExpression.Gene.10x} function.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{adhocGene.TabulateExpression.Gene.10x} function.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocGene$DE$Gene$Data}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @import Matrix
#'
#' @export
#'
#' @examples
#'
#' marvel.demo.10x <- readRDS(system.file("extdata/data",
#'                                "marvel.demo.10x.rds",
#'                                package="MARVEL")
#'                                )
#'
#' marvel.demo.10x <- adhocGene.DE.Gene.10x(MarvelObject=marvel.demo.10x)
#'
#' # Check output
#' marvel.demo.10x$adhocGene$DE$Gene$Data

adhocGene.DE.Gene.10x <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$adhocGene$Expression$Gene$Table
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$adhocGene$Expression$Gene$Table
    
    ##################################################################
        
    # Generate all possible combinations
    groups <- factor(levels(df$group), levels=levels(df$group))
    pairs <- gtools::combinations(n=length(groups), r=2, v=as.numeric(groups), repeats.allowed=FALSE)
    pairs <- as.data.frame(pairs)

    # Replace numeric values with actual factor levels
    levels <- levels(groups)

    for(i in 1:length(levels)) {
        
        pairs[pairs==i] <- levels[i]
        
    }
    
    # Compute log2FC
    log2fc <- NULL
    
    for(i in 1:nrow(pairs)) {
        
        # Retrieve cell groups
        g1 <- pairs[i,1]
        g2 <- pairs[i,2]
        
        # Retrieve mean expr
        mean.g1 <- df[which(df$group==g1), "mean.expr"]
        mean.g2 <- df[which(df$group==g2), "mean.expr"]
        
        # Compute log2FC
        log2fc[i] <- mean.g2 - mean.g1
        
    }
    
    results <- data.frame("group.pair"=paste(pairs[,2], " vs ", pairs[,1], sep=""),
                          "log2fc"=log2fc,
                          stringsAsFactors=FALSE
                          )
 
    ##################################################################
 
    # Save into new slot
    MarvelObject$adhocGene$DE$Gene$Data <- results
    
    # Return final object
    return(MarvelObject)
            
}


