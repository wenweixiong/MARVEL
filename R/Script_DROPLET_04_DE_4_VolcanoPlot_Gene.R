#' @title Plot differential gene analysis results
#'
#' @description Volcano plot of results from differential gene expression analysis. x-axis represents the log2 fold change between two cell groups. y-axis represents -log10(adjusted p-value). Only genes whose splice junctions were considered to be differentially spliced are included for plotting.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.Genes.10x} function.
#' @param pval.sj Numeric value. p-value from differential splicing analysis, below which, the splice junction is considered differentially spliced. Default is \code{0.05}.
#' @param log2fc.sj Numeric value. Absolute log2 fold change from differential splicing analysis, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{delta.sj} has been specified.
#' @param delta.sj Numeric value. Absolute difference in average PSI values between the two cell groups, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{log2fc.sj} has been specified.
#' @param min.gene.norm Numeric value. The average normalised gene expression across the two cell groups above which the splice junction is considered differentially spliced. Default is \code{0}.
#' @param pval.adj.gene Numeric value. Adjusted p-value from differential gene expression analysis, below which, the gene is considered differentially expressed. Default is \code{0.05}.
#' @param log2fc.gene Numeric value. Absolute log2 fold change from differential gene expression analysis, above which, the gene is considered differentially expressed. This option should be \code{NULL} if \code{delta.sj} has been specified.
#' @param anno Logical value. If set to \code{TRUE}, user-specific genes in \code{anno.gene_short_name} will be annotated on the plot. Default is \code{FALSE}.
#' @param anno.gene_short_name Vector of character strings. If \code{anno} set to \code{TRUE}, genes specified here will be annotated on the plot.
#' @param label.size Numeric value. If \code{anno} set to \code{TRUE}, the font size of the annotations on the plot will be adjusted to the size specified here. Default is \code{2}.
#'
#' @return An object of class S3 with a new slots \code{MarvelObject$DE$SJ$VolcanoPlot$Gene$Plot} and \code{MarvelObject$DE$SJ$VolcanoPlot$Gene$Data}.
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
#' marvel.demo.10x <- PlotDEValues.Genes.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         pval.sj=0.05,
#'                         delta.sj=5,
#'                         min.gene.norm=1.0,
#'                         pval.adj.gene=0.05,
#'                         log2fc.gene=0.5
#'                         )
#'
#' # Check outputs
#' marvel.demo.10x$DE$SJ$VolcanoPlot$Gene$Plot
#' head(marvel.demo.10x$DE$SJ$VolcanoPlot$Gene$Data)

PlotDEValues.Genes.10x <- function(MarvelObject, pval.sj=0.05, log2fc.sj=NULL, delta.sj=5, min.gene.norm=0, pval.adj.gene=0.05, log2fc.gene=0.5, anno=FALSE, anno.gene_short_name=NULL, label.size=2) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$DE$SJ$Table
    pval.sj <- pval.sj
    log2fc.sj <- log2fc.sj
    delta.sj <- delta.sj
    pval.adj.gene <- pval.adj.gene
    log2fc.gene <- log2fc.gene
    anno <- anno
    anno.gene_short_name <- anno.gene_short_name
    label.size <- label.size
    min.gene.norm <- min.gene.norm
    
    # Define arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$SJ$Table
    #pval.sj <- 0.05
    #log2fc.sj <- NULL
    #delta.sj <- 5
    #pval.adj.gene <- 0.05
    #log2fc.gene <- 0.5
    #anno <- TRUE
    #anno.gene_short_name <- df$gene_short_name[c(1:10)]
    #label.size <- 2
    #min.gene.norm <- 1
    
    ############################################################
    
    # Indicate sig events and direction
    if(!is.null(log2fc.sj)) {
        
        df$sig <- NA
        df$sig[which(df$pval < pval.sj & df$log2fc > log2fc.sj & df$mean.expr.gene.norm.g1.g2 > min.gene.norm)] <- "up"
        df$sig[which(df$pval < pval.sj & df$log2fc < (log2fc.sj*-1) & df$mean.expr.gene.norm.g1.g2 > min.gene.norm)] <- "down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
        table(df$sig)
        
    } else if(!is.null(delta.sj)){
        
        df$sig <- NA
        df$sig[which(df$pval < pval.sj & df$delta > delta.sj & df$mean.expr.gene.norm.g1.g2 > min.gene.norm)] <- "up"
        df$sig[which(df$pval < pval.sj & df$delta < (delta.sj*-1) & df$mean.expr.gene.norm.g1.g2 > min.gene.norm)] <- "down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
        table(df$sig)
        
    }
    
    # Subset sig SJ
    df <- df[which(df$sig %in% c("up", "down")), ]
    df$sig <- NULL
    
    # Keep unique gene entries
    cols.1 <- "gene_short_name"
    cols.2 <- grep("gene.norm", names(df), fixed=TRUE, value=TRUE)
    df <- df[, c(cols.1, cols.2)]
    df <- unique(df)
    
    # Report progress
    message(paste(nrow(df), " unique genes with at least one differential spliced SJ found", sep=""))
    
    
    # Indicate sig events and direction
    df$sig <- NA
    df$sig[which(df$pval.adj.gene.norm < pval.adj.gene & df$log2fc.gene.norm > log2fc.gene)] <- "up"
    df$sig[which(df$pval.adj.gene.norm < pval.adj.gene & df$log2fc.gene.norm < (log2fc.gene*-1))] <- "down"
    df$sig[is.na(df$sig)] <- "n.s."
    df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
    table(df$sig)

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
       df$label[which(df$gene_short_name %in% anno.gene_short_name)] <- df$gene_short_name[which(df$gene_short_name %in% anno.gene_short_name)]
               
    } else {
        
        df$label <- NA
        
    }
    
    ############################################################
    
    # Plot
       # Definition
       data <- df
       x <- data$log2fc.gene.norm
       y <- -log10(data$pval.adj.gene.norm)
       z <- data$sig
       label <- data$label
       maintitle <- ""
       xtitle <- "log2FC"
       ytitle <- "-log10(FDR)"
       
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
    MarvelObject$DE$SJ$VolcanoPlot$Gene$Plot <- plot
    MarvelObject$DE$SJ$VolcanoPlot$Gene$Data <- df
    
    # Return final object
    return(MarvelObject)
            
}


