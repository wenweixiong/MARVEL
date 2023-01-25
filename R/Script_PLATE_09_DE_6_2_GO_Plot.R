#' @title Plot pathway enrichment analysis results
#'
#' @description Plots user-specified enriched pathways.
#'
#' @details
#' This function plots selected gene sets returned from gene ontology analysis performed previously using \code{BioPathways}
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{BioPathways} function.
#' @param go.terms Vector of character strings. Names of pathways to plot. Should match pathway names in column \code{Description} of \code{MarvelObject$DE$BioPathways$Table}.
#' @param y.label.size Numeric value. Size of y-axis tick labels, i.e. gene set names.
#' @param offset Numeric value. The -log10(p-value) on the x-axis to substract or add to increase the plot margins.
#' @param x.axis Character string. If set to \code{"enrichment"} (default) the pathway enrichment will be displayed on the x-axis while the color intensity of the data points will reflect the -log10(adjusted p-value). If set to \code{"pval"} the -log10(adjusted p-value) will be displayed on the x-axis while the color intensity of the data points will reflect the pathway enrichment.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$DE$BioPathways$Plot}.
#'
#' @importFrom plyr join
#' @import methods
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' # Define go terms to plot
#' df <- marvel.demo$DE$BioPathways$Table
#' go.terms <- df$Description[c(1:10)]
#'
#' # Plot
#' marvel.demo <- BioPathways.Plot(MarvelObject=marvel.demo,
#'                                 go.terms=go.terms,
#'                                 offset=10
#'                                 )
#'
#' # Check output
#' marvel.demo$DE$BioPathways$Plot

BioPathways.Plot <- function(MarvelObject, go.terms, y.label.size=10, offset=0.5, x.axis="enrichment") {
    
    # Define arguments
    df <- MarvelObject$DE$BioPathways$Table
    go.terms <- go.terms
    y.label.size <- y.label.size
    offset <- offset
    x.axis <- x.axis
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- marvel$DE$BioPathways$Table
    #go.terms <- go.terms
    #y.label.size <- 10
    #offset <- 0.5
    #x.axis <- "pval"
    
    ###############################################################

    # Subset GO terms to plot
    df <- df[which(df$Description %in% go.terms), ]
    
    # transform adjusted p-val
    df$p.adjust <- -log10(df$p.adjust)

    # Order by pval/enrichment
    if(x.axis=="pval") {
        
        df <- df[order(df$p.adjust, decreasing=TRUE), ]
        
    } else if(x.axis=="enrichment"){
        
        df <- df[order(df$enrichment, decreasing=TRUE), ]
        
    }
    
    # Wrap label
    labels <- data.frame("x"=as.character(df$Description), stringsAsFactors=FALSE)
    labels$y <- ifelse(nchar(labels$x) < 30, 10, 20)
    labels$newx <- stringr::str_wrap(labels$x, width=25)
    df$Description <- labels$newx
    df$Description <- factor(df$Description, levels=rev(df$Description))
    
    # Dotplot
    if(x.axis=="pval") {
        
        # Definitions
        data <- df
        x <- data$Description
        z2 <- data$enrichment
        y <- data$p.adjust
        z1 <- data$Count
        maintitle <- ""
        xtitle <- ""
        legendtitle.color <- "Pathway enrichment"
        #fivenum(y) ; ymin <- 3.5 ; ymax <- 5.5 ; yinterval <- 0.5
        ytitle  <- "-log10(p-value)"
        legendtitle.size <- "Gene hits"
        
        ymin <- min(y) - offset ; ymax <- max(y) + offset
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, color=z2, size=z1)) +
            scale_y_continuous(limits=c(ymin, ymax)) +
            scale_color_gradient(low="blue", high="red") +
            labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle.color, size=legendtitle.size) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black"),
                axis.text.y=element_text(size=y.label.size, colour="black"),
                legend.title=element_text(size=10),
                legend.text=element_text(size=10)
                ) +
            coord_flip()
        
        
    } else if(x.axis=="enrichment"){
        
        # Definitions
        data <- df
        x <- data$Description
        y <- data$enrichment
        z2 <- data$p.adjust
        z1 <- data$Count
        maintitle <- ""
        xtitle <- ""
        ytitle <- "Pathway enrichment"
        #fivenum(y) ; ymin <- 3.5 ; ymax <- 5.5 ; yinterval <- 0.5
        legendtitle.color <- "-log10(p-value)"
        legendtitle.size <- "Gene hits"
        
        ymin <- min(y) - offset ; ymax <- max(y) + offset
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, color=z2, size=z1)) +
            scale_y_continuous(limits=c(ymin, ymax)) +
            scale_color_gradient(low="blue", high="red") +
            labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle.color, size=legendtitle.size) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black"),
                axis.text.y=element_text(size=y.label.size, colour="black"),
                legend.title=element_text(size=10),
                legend.text=element_text(size=10)
                ) +
            coord_flip()
            
    }
                
    # Save into new slot
    MarvelObject$DE$BioPathways$Plot <- plot
       
    return(MarvelObject)
    
}
