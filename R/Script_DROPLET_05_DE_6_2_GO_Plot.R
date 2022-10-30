#' @title Plot pathway enrichment analysis results
#'
#' @description Plots user-specified enriched pathways.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{BioPathways.10x} function.
#' @param go.terms Vector of character strings. Names of pathways to plot. Should match pathway names in column \code{Description} of \code{MarvelObject$DE$BioPathways$Table}.
#' @param y.label.size Numeric value. Size of y-axis tick labels, i.e. pathway names.
#' @param offset Numeric value. The value on the x-axis to substract or add to increase the plot margins.
#'
#' @return An object of class S3 containing with new slot \code{MarvelObject$DE$BioPathways$Plot}.
#'
#' @importFrom plyr join
#' @import methods
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
#' # Define top pathways to plot
#' go.terms <- marvel.demo.10x$DE$BioPathways$Table$Description
#' go.terms <- go.terms[c(1:10)]
#'
#' # Plot
#' marvel.demo.10x <- BioPathways.Plot.10x(
#'                             MarvelObject=marvel.demo.10x,
#'                             go.terms=go.terms
#'                             )
#' # Check outpout
#' marvel.demo.10x$DE$BioPathways$Plot

BioPathways.Plot.10x <- function(MarvelObject, go.terms, y.label.size=10, offset=0.5) {
    
    # Define arguments
    df <- MarvelObject$DE$BioPathways$Table
    go.terms <- go.terms
    y.label.size <- y.label.size
    offset <- offset
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$BioPathways$Table
    #go.terms <- go.terms
    #y.label.size <- 10
    #offset <- 0.5
    
    # Subset GO terms to plot
    df <- df[which(df$Description %in% go.terms), ]
    
    # Transform adjusted p-val
    df$p.adjust <- -log10(df$p.adjust)
    
    # Order by enrichment
    df <- df[order(df$enrichment, decreasing=TRUE), ]

    # Wrap label
    labels <- data.frame("x"=as.character(df$Description), stringsAsFactors=FALSE)
    labels$y <- ifelse(nchar(labels$x) < 30, 10, 20)
    labels$newx <- stringr::str_wrap(labels$x, width=25)
    df$Description <- labels$newx
    df$Description <- factor(df$Description, levels=rev(df$Description))
    
    # Dotplot
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
                
    # Save into new slot
    MarvelObject$DE$BioPathways$Plot <- plot
       
    return(MarvelObject)
    
}
