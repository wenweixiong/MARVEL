#' @title Plot differential splice junction analysis results
#'
#' @description Volcano plot of results from differential splice junction analysis. x-axis represents the average normalised gene expression across the two cell groups. y-axis represents the differences or log2 fold change between the two cell groups.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.Genes.10x} function.
#' @param pval Numeric value. p-value, below which, the splice junction is considered differentially spliced. To be used in conjunction with \code{log2fc}, \code{delta}, and \code{min.gene.norm}. Default is \code{0.05}.
#' @param log2fc Numeric value. Absolute log2 fold change, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{delta} has been specified.
#' @param delta Numeric value. Absolute differences in average PSI values between the two cell groups, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{log2fc} has been specified.
#' @param min.gene.norm Numeric value. The average normalised gene expression across the two cell groups above which the splice junction is considered differentially spliced. Default is \code{0}.
#' @param anno Logical value. If set to \code{TRUE}, user-specific spliced genes in \code{anno.coord.intron} will be annotated on the plot. Default is \code{FALSE}.
#' @param anno.coord.intron Vector of character strings. If \code{anno} set to \code{TRUE}, splice junction coordinates specified here will be annotated on the plot.
#' @param label.size Numeric value. If \code{anno} set to \code{TRUE}, the font size of the annotations on the plot will be adjusted to the size specified here. Default is \code{2}.
#'
#' @return An object of class S3 with a new slots \code{MarvelObject$DE$SJ$VolcanoPlot$SJ$Plot} and \code{MarvelObject$DE$SJ$VolcanoPlot$SJ$Data}.
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
#' marvel.demo.10x <- PlotDEValues.SJ.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         pval=0.05,
#'                         delta=5,
#'                         min.gene.norm=1.0,
#'                         anno=FALSE
#'                         )
#'
#' # Check outputs
#' marvel.demo.10x$DE$SJ$VolcanoPlot$SJ$Plot
#' head(marvel.demo.10x$DE$SJ$VolcanoPlot$SJ$Data)

PlotDEValues.SJ.10x <- function(MarvelObject, pval=0.05, log2fc=NULL, delta=5, min.gene.norm=0, anno=FALSE, anno.coord.intron=NULL, label.size=2) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$DE$SJ$Table
    pval <- pval
    log2fc <- log2fc
    delta <- delta
    anno <- anno
    anno.coord.intron <- anno.coord.intron
    label.size <- label.size
    min.gene.norm <- min.gene.norm
    
    # Define arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$SJ$Table
    #pval <- 0.05
    #log2fc <- NULL
    #delta <- 5
    #anno <- FALSE
    #anno.coord.intron <- df$coord.intron[c(1:10)]
    #label.size <- 2
    #min.gene.norm <- 1
    
    ############################################################
    
    # Subset/plot expressed genes
    df <- df[which(df$mean.expr.gene.norm.g1.g2 > min.gene.norm), ]
    
    # Indicate sig events and direction
    if(!is.null(log2fc)) {
        
        df$sig <- NA
        df$sig[which(df$pval < pval & df$log2fc > log2fc)] <- "up"
        df$sig[which(df$pval < pval & df$log2fc < (log2fc*-1))] <- "down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
        table(df$sig)
        
    } else if(!is.null(delta)){
        
        df$sig <- NA
        df$sig[which(df$pval < pval & df$delta > delta)] <- "up"
        df$sig[which(df$pval < pval & df$delta < (delta*-1))] <- "down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
        table(df$sig)
        
    }
    
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
       df$label[which(df$coord.intron %in% anno.coord.intron)] <- df$gene_short_name[which(df$coord.intron %in% anno.coord.intron)]
               
    } else {
        
        df$label <- NA
        
    }

    # Compute mean norm gene expression
    numerator <- (df$mean.expr.gene.norm.g1 * df$n.cells.expr.gene.norm.g1) + (df$mean.expr.gene.norm.g2 * df$n.cells.expr.gene.norm.g2)
    denominator <- df$n.cells.expr.gene.norm.g1 + df$n.cells.expr.gene.norm.g2
    df$mean.expr.gene.norm.g1.g2 <- numerator/denominator
    
    ############################################################
    
    # Plot
       # Definition
       data <- df
       x <- data$mean.expr.gene.norm.g1.g2
       y <- data$delta
       z <- data$sig
       label <- data$label
       maintitle <- ""
       xtitle <- "Mean log2(Norm. GEX + 1)"
       
       if(!is.null(log2fc)) {
           
           ytitle <- "log2fc"
           
       } else if(!is.null(delta)){
           
           ytitle <- "delta PSI"
           
       }
            
       # Plot
       plot <- ggplot() +
                  geom_point(data, mapping=aes(x=x, y=y, color=z), shape=20, alpha = 0.75, size=0.1) +
                  ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                  scale_colour_manual(values=col.breaks) +
                  #scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                  #scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
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
                           
    # Save to new slot
    MarvelObject$DE$SJ$VolcanoPlot$SJ$Plot <- plot
    MarvelObject$DE$SJ$VolcanoPlot$SJ$Data <- df
    
    # Return final object
    return(MarvelObject)
            
}


