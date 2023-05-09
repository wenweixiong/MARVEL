#' @title Elbow plot
#'
#' @description Create elbw plot from principal components (PCs).
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{RunPCA} function.
#' @param level Character strings. Indicate \code{"splicing"} or \code{"gene"} to retrieve the eigenvalues from splicing or gene expression principal component analysis (PCA), respectively.
#' @param n.dim. Indicate the first number of dimensions to plot the eigenvalues for. Default value is \code{50}, i.e., the first 50 PCs.
#'
#' @return An object of class S3 containing with new slots \code{ MarvelObject$PCA$PSI$ElbowPlot} or \code{MarvelObject$PCA$Exp$ElbowPlot} when \code{level} option set to \code{"splicing"} or \code{"gene"}, respectively.
#'
#' @importFrom plyr join
#' @import methods
#' @import ggplot2
#'
#' @export
#'

ElbowPlot <- function(MarvelObject, level, n.dim=50) {

    # Define arguments
    MarvelObject <- MarvelObject
    level <- level
    n.dim <- n.dim
    
    # Example arguments
    #MarvelObject <- marvel
    #level <- "gene"
    #n.dim <- 50
    
    ######################################################################
    
    # Retrieve eigen values
    if(level=="splicing") {
       
       df <- MarvelObject$PCA$PSI$EigenValue
        
    } else if(level=="gene"){
        
        df <- MarvelObject$PCA$Exp$EigenValue
        
    }
    
    # Subset dimensions
    index <- c(1:n.dim)
    df <- df[index, ]
     
    # Dot plot
        # Definition
        data <- df
        x <- data$pc
        y <- data$variance.percent
        maintitle <- ""
        xtitle <- "PC"
        ytitle <- "Variance explained (%)"
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y), fill="black", size=1) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                plot.title = element_text(size=12, hjust=0.5),
                axis.line=element_line(colour = "black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black"),
                axis.text.y=element_text(size=10, colour="black")
                )

    ######################################################################
    
    # Save to new slot
    if(level=="splicing") {
        
        MarvelObject$PCA$PSI$ElbowPlot <- plot

    } else if(level=="gene") {
        
        MarvelObject$PCA$Exp$ElbowPlot <- plot
        
    }

    return(MarvelObject)
        
}
