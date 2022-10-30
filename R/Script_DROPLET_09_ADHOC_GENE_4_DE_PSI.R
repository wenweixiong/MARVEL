#' @title Differential splice junction analysis of specified gene
#'
#' @description Performs differential splice junction analysis specified gene across for all possible pairs of cell groups. The gene and cell groups were defined earlier in \code{adhocGene.TabulateExpression.Gene.10x} function.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{adhocGene.TabulateExpression.PSI.10x} function.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocGene$DE$PSI$Data}.
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
#' marvel.demo.10x <- adhocGene.DE.PSI.10x(MarvelObject=marvel.demo.10x)
#'
#' # Check output
#' marvel.demo.10x$adhocGene$DE$PSI$Data

adhocGene.DE.PSI.10x <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$adhocGene$Expression$PSI$Table
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$adhocGene$Expression$PSI$Table
    
    ########################################################################
    
    # Generate all possible combinations
    groups <- factor(levels(df$group), levels=levels(df$group))
    pairs <- gtools::combinations(n=length(groups), r=2, v=as.numeric(groups), repeats.allowed=FALSE)
    pairs <- as.data.frame(pairs)

    # Replace numeric values with actual factor levels
    levels <- levels(groups)

    for(i in 1:length(levels)) {
        
        pairs[pairs==i] <- levels[i]
        
    }
    
    # Compute delta PSI
    coord.introns <- levels(df$coord.intron)
    
    .list <- list()
    
    for(j in 1:length(coord.introns)) {
        
        df.small <- df[which(df$coord.intron==coord.introns[j]), ]
        
        if(nrow(df.small) >=2) {
            
            pairs.small <- pairs[which(pairs$V2 %in% df.small$group), ]
            pairs.small <- pairs.small[which(pairs.small$V1 %in% df.small$group), ]
            
            delta <- NULL
            
            for(i in 1:nrow(pairs.small)) {
                
                # Retrieve cell groups
                g1 <- pairs.small[i,1]
                g2 <- pairs.small[i,2]
                
                # Retrieve mean expr
                mean.g1 <- df.small[which(df.small$group==g1), "psi"]
                mean.g2 <- df.small[which(df.small$group==g2), "psi"]
                
                # Compute log2FC
                delta[i] <- mean.g2 - mean.g1
                
            }
            
            results <- data.frame("coord.intron"=coord.introns[j],
                                  "group.pair"=paste(pairs.small[,2], " vs ", pairs.small[,1], sep=""),
                                  "delta"=delta,
                                  stringsAsFactors=FALSE
                                  )
                                  
            .list[[j]] <- results
                        
        }
        
    }
    
    results <- do.call(rbind.data.frame, .list)

    # Annotate column no.
        # Create reference table
        . <- unique(df[,c("coord.intron", "figure.column")])
        
        # Annotate
        results <- join(results, ., by="coord.intron", type="left")
        
        # Reorder columns
        cols <- c("group.pair", "figure.column", "coord.intron", "delta")
        results <- results[, cols]
        
    ########################################################################
    
    # Save into new slots
    MarvelObject$adhocGene$DE$PSI$Data <- results
    
    # Return final object
    return(MarvelObject)
            
}


