#' @title Plot differential splicing analysis results based on mean PSI difference
#'
#' @description Volcano plot of differential splicing analysis results based on mean PSI difference between 2 groups of cells. x-axis represents the mean delta PSI. y-axis represents the adjusted p-values.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues} function.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param pval Numeric value. Adjusted p-value below which the splcing event are considered as statistically significant and will consequently be color-annotated on the plot.
#' @param delta Numeric value. The positive (and negative) value specified above (and below) which the splicing events are considered to be statistically significant and will consequently be color-annotated on the plot.
#' @param point.size Numeric value. The point size for the data points. Default value is \code{1}.
#' @param anno Logical value. If set to \code{TRUE}, the specific gene names will be annotated on the plot. Speficified together with \code{anno.tran_id}.
#' @param anno.tran_id Vector of character strings. When \code{anno} set to \code{TRUE}, the coordinates of the splicing events to be annotated on the plot.
#' @param label.size Numeric value. Only applicable if \code{anno} set to TRUE. Size of the gene name labels.
#' @param y.upper.offset Numeric value. The value in -log10(p-value) to increase the upper limit of the y-axis. To be used when \code{anno} set to TRUE so that gene labels will not be truncated at the upper limit of the y-axis.
#' @param xlabel.size Numeric value. Font size of the xtick labels. Default is \code{8}.
#'
#' @return An object of class S3 containing with new slot \code{MarvelObject$DE$PSI$Plot[["method"]]}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @import scales
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- PlotDEValues.PSI.Mean(MarvelObject=marvel.demo,
#'                                      method="ad",
#'                                      pval=0.10,
#'                                      delta=5
#'                                      )
#'
#' # Check output
#' marvel.demo$DE$PSI$Plot[["ad"]]

PlotDEValues.PSI.Mean <- function(MarvelObject, method, pval=0.1, delta=5, point.size=1, xlabel.size=8, anno=FALSE, anno.tran_id=NULL, label.size=2.5, y.upper.offset=5) {

    # Define arguments
    MarvelObject <- MarvelObject
    method <- method
    df <- MarvelObject$DE$PSI$Table[[method]]
    pval <- pval
    delta <- delta
    point.size <- point.size
    anno <- anno
    anno.tran_id <- anno.tran_id
    label.size <- label.size
    y.upper.offset <- y.upper.offset
    xlabel.size <- xlabel.size
    
    # Example arguments
    #MarvelObject <- marvel
    #method <- "ad"
    #df <- MarvelObject$DE$PSI$Table[[method]]
    #pval <- 0.10
    #delta <- 5
    #anno <- TRUE
    #anno.tran_id <- tran_ids
    #label.size <- 2.5
    #y.upper.offset <- 5
    #point.size <- 1
    #label.size <- 2.5
    #xlabel.size <- 8.0
    
    ##############################################
    
    # Remove outliers
    df <- df[which(df$outlier==FALSE), ]
    
    # Indicate sig events and direction
    df$sig <- NA
    df$sig[which(df$p.val.adj < pval & df$mean.diff > delta)] <- "up"
    df$sig[which(df$p.val.adj < pval & df$mean.diff < (delta*-1))] <- "down"
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
        
        df$label <- NA
        df$gene_short_name.event_type <- paste(df$gene_short_name, " (", df$event_type, ")", sep="")
                    
        index <- which(df$tran_id %in% anno.tran_id)
        df$label[index] <- df$gene_short_name.event_type[index]
            
       # Definition
       data <- df
       x <- data$mean.diff
       y <- -log10(data$p.val.adj)
       z <- data$sig
       label <- data$label
       maintitle <- ""
       xtitle <- "delta PSI"
       ytitle <- "-log10(p-value)"

       xmin <- -100 ; xmax <- 100 ; xinterval <- 25
       ymin <- 0 ; ymax <- max(y) + y.upper.offset ; yinterval <- 5
       
       # Plot
       plot <- ggplot() +
                  geom_point(data, mapping=aes(x=x, y=y, color=z), shape=20, alpha = 0.75, size=point.size) +
                  ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
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
                       axis.text.x=element_text(size=xlabel.size, colour="black"),
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

        xmin <- -100 ; xmax <- 100 ; xinterval <- 25
        ymin <- 0 ; ymax <- max(y) ; yinterval <- 5
 
        # Plot
        plot <- ggplot() +
                   geom_point(data, mapping=aes(x=x, y=y, color=z), shape=20, alpha = 0.75, size=point.size) +
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
                        axis.text.x=element_text(size=xlabel.size, colour="black"),
                        axis.text.y=element_text(size=10, colour="black"),
                        legend.position="none",
                        legend.title=element_text(size=8),
                        legend.text=element_text(size=8)
                        )
    
    }
                 
    ##############################################
    
    # Save to new slot
    MarvelObject$DE$PSI$Plot[[method]] <- plot
  
    # Return final object
    return(MarvelObject)
        
}
