#' @title Plot global differential gene expression analysis results
#'
#' @description Volcano plot of differential splicing analysis results based on all expressed genes between 2 groups of cells. x-axis represents the log2 fold change in gene expression. y-axis represents the adjusted p-values.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues} function.
#' @param pval Numeric value. Adjusted p-value below which the genes are considered as statistically significant and will consequently be color-annotated on the plot.
#' @param log2fc Numeric value. The positive (and negative) value specified above (and below) which the genes are considered to be statistically significant and will consequently be color-annotated on the plot.
#' @param point.size Numeric value. The point size for the data points. Default value is \code{1}.
#' @param anno Logical value. If set to \code{TRUE}, the specific gene names will be annotated on the plot as defined in \code{anno.gene_short_name} option.
#' @param anno.gene_short_name Vector of character strings. When \code{anno} set to \code{TRUE}, the gene names to be annotated on the plot.
#' @param label.size Numeric value. Only applicable if \code{anno} set to TRUE. Size of the gene name labels.
#' @param y.upper.offset Numeric value. The value in -log10(p-value) to increase the upper limit of the y-axis. To be used when \code{anno} set to TRUE so that gene labels will not be truncated at the upper limit of the y-axis.
#' @param xlabel.size Numeric value. Font size of the xtick labels. Default is \code{8}.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$DE$Exp.Global$Table}, \code{MarvelObject$DE$Exp.Global$Summary}, and \code{MarvelObject$DE$Exp.Global$Plot}
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
#' marvel.demo <- PlotDEValues.Exp.Global(MarvelObject=marvel.demo,
#'                                        pval=0.10,
#'                                        log2fc=0.5
#'                                        )
#'
#' # Check output
#' head(marvel.demo$DE$Exp.Global$Table)
#' marvel.demo$DE$Exp.Global$Plot
#' marvel.demo$DE$Exp.Global$Summary

PlotDEValues.Exp.Global <- function(MarvelObject, pval=0.10, log2fc=0.5, point.size=1, anno=FALSE, anno.gene_short_name=NULL, label.size=2.5, y.upper.offset=5, xlabel.size=8) {

    # Define arguments
    df <- MarvelObject$DE$Exp$Table
    pval <- pval
    log2fc <- log2fc
    point.size <- point.size
    anno <- anno
    anno.gene_short_name <- anno.gene_short_name
    y.upper.offset <- y.upper.offset
    xlabel.size <- xlabel.size
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$Exp$Table
    #pval <- 0.10
    #log2fc <- 0.5
    #anno <- FALSE
    #anno.gene_short_name <- gene_short_names
    #label.size <- 2.5
    #point.size <- 1
    #y.upper.offset <- 5
    #xlabel.size <- 2.5
    
    # Indicate sig events and direction
    df$sig <- NA
    df$sig[which(df$p.val.adj < pval & df$log2fc > log2fc)] <- "up"
    df$sig[which(df$p.val.adj < pval & df$log2fc < (log2fc*-1))] <- "down"
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
        
       df$label <- ifelse(df$gene_short_name %in% anno.gene_short_name, df$gene_short_name, "")
                           
       # Definition
       data <- df
       x <- data$log2fc
       y <- -log10(data$p.val.adj)
       z <- data$sig
       label <- data$label
       maintitle <- ""
       xtitle <- "log2fc"
       ytitle <- "-log10(p-value)"

       xmin <- floor(min(x)) ; xmax <- ceiling(max(x)) ; xinterval <- 2
       ymin <- 0 ; ymax <- max(y) + y.upper.offset ; yinterval <- 5
       
       if((xmin %% 2) != 0) {
           
           xmin <- xmin - 1 ; xmax <- xmax + 1
           
       }
       
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
        x <- data$log2fc
        y <- -log10(data$p.val.adj)
        z <- data$sig
        maintitle <- ""
        xtitle <- "log2FC"
        ytitle <- "-log10(p-value)"

        xmin <- floor(min(x)) ; xmax <- ceiling(max(x)) ; xinterval <- 2
        ymin <- 0 ; ymax <- max(y) ; yinterval <- 5
 
        if((xmin %% 2) != 0) {
             
             xmin <- xmin - 1 ; xmax <- xmax + 1
             
        }
 
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
    
    # Summary
    tbl <- as.data.frame(table(df$sig))
    names(tbl) <- c("sig", "freq")
    
    ################################################

    # Save to new slot
    MarvelObject$DE$Exp.Global$Summary <- tbl
    MarvelObject$DE$Exp.Global$Table <- df
    MarvelObject$DE$Exp.Global$Plot <- plot

    # Return final object
    return(MarvelObject)
        
}
