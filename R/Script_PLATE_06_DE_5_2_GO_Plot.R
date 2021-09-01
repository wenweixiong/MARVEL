#' @title Plot Gene Ontology Analysis Results
#'
#' @description
#' \code{BioPathways.Plot} plots selected gene sets returned from gene ontology analysis.
#'
#' @details
#' This function plots selected gene sets returned from gene ontology analysis performed previously using \code{BioPathways}
#'
#' @param MarvelObject S3 object generated from \code{CompareValues} function.
#' @param go.terms Character strings. Names of gene sets to plot. Should match gene sets name in column \code{Term} of \code{MarvelObject$DE$BioPathways}.
#' @param y.label.size Numeric value. Size of y-axis tick labels, i.e. gene set names.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to a new slot named \code{$DE$BioPathways.Plot} containing a plot of gene set name on the y-axis, odds ratio on the x-axis, data point color gradient by adjusted p-values, and data point size scaled by number of differentially spliced genes found in each gene set.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @import stringr
#' @importFrom AnnotationDbi select
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Specify GO terms to plot
#' go.terms <- marvel$DE$BioPathways$Table$Term[c(1:10)]
#'
#' # Run example
#' marvel <- BioPathways.Plot(MarvelObject=marvel,
#'                            go.terms=go.terms,
#'                            y.label.size=10
#'                            )
#'
#' # Check output
#' marvel$DE$BioPathways$Plot

BioPathways.Plot <- function(MarvelObject, go.terms, y.label.size=10) {
    
    # Define arguments
    df <- MarvelObject$DE$BioPathways$Table
    go.terms <- go.terms
    y.label.size <- y.label.size
    
    # Example arguments
    #df <- marvel$DE$BioPathways$Table
    #go.terms <-
    #y.label.size <- 10
    
    # Re-code Inf OR
        # Retrieve largest non-Inf value
        #. <- sort(df$OddsRatio[c(1:10)][which(df$OddsRatio[c(1:10)] != Inf)], decreasing=TRUE)[1]
        
        # Re-code Inf
        df$OddsRatio[which(df$OddsRatio==Inf)] <- 10
    
    # Subset GO terms to plot
    df <- df[which(df$Term %in% go.terms), ]
    
    # transform adjusted p-val
    df$p.val.adj <- -log10(df$p.val.adj)
    
    # Order by odds ratio
    df <- df[order(df$OddsRatio, decreasing=TRUE), ]

    # Wrap label
    labels <- data.frame("x"=as.character(df$Term), stringsAsFactors=FALSE)
    labels$y <- ifelse(nchar(labels$x) < 25, 10, 20)
    labels$newx <- str_wrap(labels$x, width=25)
    df$Term <- labels$newx
    df$Term <- factor(df$Term, levels=rev(df$Term))
    
    # Dotplot
        # Definitions
        data <- df
        x <- data$Term
        y <- data$OddsRatio
        z1 <- data$p.val.adj
        z2 <- data$Count
        maintitle <- ""
        xtitle <- ""
        ytitle <- "Enrichment (Gene Ratio)"
        #fivenum(y) ; ymin <- 3.5 ; ymax <- 5.5 ; yinterval <- 0.5
        legendtitle.color <- "-log10(p-val)"
        legendtitle.size <- "Counts"
        
        ymin <- min(y) - 1 ; ymax <- max(y)
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, color=z1, size=z2)) +
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
