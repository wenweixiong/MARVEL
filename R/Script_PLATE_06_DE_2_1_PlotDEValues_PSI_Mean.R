#' @title Plot Differential Splicing Analysis Results
#'
#' @description
#' \code{PlotDEValues.PSI.Mean} plots the differential splicing analysis results.
#'
#' @details
#' This function plot the differential splicing analysis results and can annotate the user-defined significant data points and gene names. Typically used when the means of percent spliced-in (PSI) values of 2 groups of cells were previously tested using t-test or Wilcoxon rank-sum test.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues} function.
#' @param de.p.val.adj Numeric value. Adjusted p-value below which the splcing event are considered as statistically significant and will consequently be color-annotated on the plot.
#' @param de.mean.diff Numeric value. The positive (and negative) value specified above (and below) which the splicing events are considered to be statistically significant and will consequently be color-annotated on the plot.
#' @param anno Logical value. If set to \code{TRUE}, the specific gene names will be annotated on the plot. Speficified together with \code{anno.x.pos}, \code{anno.y.pos}, \code{anno.x.neg}, and \code{anno.y.neg}.
#' @param anno.x.pos Numeric value. The value above on the x-axis which the gene names will be annotated on the plot. Specified together with \code{anno.y.pos}, \code{anno.x.neg}, and \code{anno.y.neg}.
#' @param anno.y.pos Numeric value. The value above on the y-axis which the gene names will be annotated on the plot. Specified together with \code{anno.x.pos}, \code{anno.x.neg}, and \code{anno.y.neg}.
#' @param anno.x.neg Numeric value. The value below on the x-axis which the gene names will be annotated on the plot. Specified together with \code{anno.x.pos}, \code{anno.y.pos}, and \code{anno.y.neg}.
#' @param anno.y.neg Numeric value. The value above on the y-axis which the gene names will be annotated on the plot. Specified together with \code{anno.y.pos}, \code{anno.x.neg}, and \code{anno.x.pos}.
#' @param label.size Numeric value. Only applicable if \code{anno} set to TRUE. Size of the gene name labels.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$DE$PSI$Plot}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import ggplot2
#' @import ggrepel
#' @import scales
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- PlotDEValues.PSI.Mean(MarvelObject=marvel,
#'                                 de.p.val.adj=0.10,
#'                                 de.mean.diff=0.01,
#'                                 anno=FALSE
#'                                 )
#'
#' # Check output
#' marvel$DE$PSI$Plot

PlotDEValues.PSI.Mean <- function(MarvelObject, de.p.val.adj, de.mean.diff, anno=FALSE, anno.x.pos=NULL, anno.y.pos=NULL, anno.x.neg=NULL, anno.y.neg=NULL, label.size=2.5) {

    # Define arguments
    df <- MarvelObject$DE$PSI$Table
    de.p.val.adj <- de.p.val.adj
    de.mean.diff <- de.mean.diff
    anno <- anno
    anno.x.pos <- anno.x.pos
    anno.y.pos <- anno.y.pos
    anno.x.neg <- anno.x.neg
    anno.y.neg <- anno.y.neg
    label.size <- label.size
    
    # Example arguments
    #df <- marvel$DE$PSI
    #de.p.val.adj <- 0.10
    #de.mean.diff <- 0.05
    #anno <- TRUE
    #anno.x.pos <- 0.10
    #anno.y.pos <- 15
    #anno.x.neg <- -0.10
    #anno.y.neg <- 15
    #label.size <- 2.5
    
    # Indicate sig events and direction
    df$sig <- NA
    df$sig[which(df$p.val.adj < de.p.val.adj & df$mean.diff > de.mean.diff)] <- "up"
    df$sig[which(df$p.val.adj < de.p.val.adj & df$mean.diff < (de.mean.diff*-1))] <- "down"
    df$sig[is.na(df$sig)] <- "n.s."
    df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
    
    # Indicate color scheme
    sig.up <- which(df$sig=="up")
    sig.down <- which(df$sig=="down")

    if(length(sig.up) != 0 & length(sig.down) != 0) {

        col.breaks <- c("red", "blue", "gray")
        
    } else if(length(sig.up) != 0 & length(sig.down) == 0) {

        col.breaks <- c("red", "gray")
        
    } else if(length(sig.up) == 0 & length(sig.down) != 0) {

        col.breaks <- c("blue", "gray")

    } else if(length(sig.up) == 0 & length(sig.down) == 0) {

        col.breaks <- "gray"
        
    }
    
    # Create labels
    if(anno==TRUE) {
        
        df$gene_short_name.event_type <- paste(df$gene_short_name, ".", df$event_type, sep="")
        df$label <- NA
        df$label[which(-log10(df$p.val.adj) > anno.y.pos & df$mean.diff > anno.x.pos)] <- df$gene_short_name.event_type[which(-log10(df$p.val.adj) > anno.y.pos & df$mean.diff > anno.x.pos)]
        df$label[which(-log10(df$p.val.adj) > anno.y.neg & df$mean.diff < anno.x.neg)] <- df$gene_short_name.event_type[which(-log10(df$p.val.adj) > anno.y.neg & df$mean.diff < anno.x.neg)]
                           
       # Definition
       data <- df
       x <- data$mean.diff
       y <- -log10(data$p.val.adj)
       z <- data$sig
       label <- data$label
       maintitle <- ""
       xtitle <- "delta PSI"
       ytitle <- "-log10(p-value)"

       xmin <- -1 ; xmax <- 1 ; xinterval <- 0.25
       ymin <- 0 ; ymax <- max(y) + 5 ; yinterval <- 5
       
       # Plot
       plot <- ggplot() +
                  geom_point(data, mapping=aes(x=x, y=y, color=z), shape=20, alpha = 0.75, size=0.5) +
                  geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                  scale_colour_manual(values=col.breaks) +
                  scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                  scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                  labs(title=maintitle, x=xtitle, y=ytitle) +
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
                       axis.text.y=element_text(size=10, colour="black"),
                       legend.position="none",
                       legend.title=element_text(size=8),
                       legend.text=element_text(size=8)
                       )
    
    } else {
        
        # Definition
        data <- df
        x <- data$mean.diff
        y <- -log10(data$p.val.adj)
        z <- data$sig
        maintitle <- ""
        xtitle <- "delta PSI"
        ytitle <- "-log10(p-value)"

        xmin <- -1 ; xmax <- 1 ; xinterval <- 0.25
        ymin <- 0 ; ymax <- max(y) ; yinterval <- 5
 
        # Plot
        plot <- ggplot() +
                   geom_point(data, mapping=aes(x=x, y=y, color=z), shape=20, alpha = 0.75, size=0.5) +
                   scale_colour_manual(values=col.breaks) +
                   scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                   scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                   labs(title=maintitle, x=xtitle, y=ytitle) +
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
                        axis.text.y=element_text(size=10, colour="black"),
                        legend.position="none",
                        legend.title=element_text(size=8),
                        legend.text=element_text(size=8)
                        )
    
    }
                          
    # Save to new slot
    MarvelObject$DE$PSI$Plot <- plot
  
    return(MarvelObject)
        
}
