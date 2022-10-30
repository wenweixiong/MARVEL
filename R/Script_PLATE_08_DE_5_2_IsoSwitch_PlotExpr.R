#' @title Plot gene-splicing relative change
#'
#' @description Plots delta PSI vs gene log2-fold change
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.Genes.10x} function.
#' @param anno Logical value. If set to \code{TRUE}, genes with coordinated, opposing or complex change relative to splicing change will be annotated on the plot. Default value is \code{FALSE}.
#'
#' @return An object of class S3 containing with new slots \code{MarvelObject$DE$Cor$PSIvsExpr$Plot}.
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
#' marvel.demo <- IsoSwitch.PlotExpr(MarvelObject=marvel.demo, anno=TRUE)
#'
#' # Check output
#' marvel.demo$DE$Cor$PSIvsExpr$Plot

IsoSwitch.PlotExpr <- function(MarvelObject, anno=FALSE) {

    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$DE$Cor$Table_Raw
    psi.delta <- MarvelObject$DE$Cor$param$psi.delta
    gene.log2fc <- MarvelObject$DE$Cor$param$gene.log2fc
    anno <- anno
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$Cor$Table_Raw
    #psi.delta <- MarvelObject$DE$Cor$param$psi.delta
    #gene.log2fc <- MarvelObject$DE$Cor$param$gene.log2fc
    #anno <- TRUE
    
    ##########################################
    
    # Set factor levels
    levels.complete <- c("Coordinated",
                         "Opposing",
                         "Iso-Switch",
                         "Complex"
                         )
    levels <- intersect(levels.complete, unique(df$cor))
    df$cor <- factor(df$cor, levels=levels)
     
    # Define color scheme
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    cols.complete = gg_color_hue(4)
    
    color.df <- data.frame("level"=levels.complete,
                           "color"=cols.complete
                           )
    color.df <- color.df[which(color.df$level %in% levels), ]
    colors <- color.df$color
    
    if(anno==FALSE) {
        
        # Definition
        data <- df
        x <- data$mean.diff
        y <- data$log2fc.gene
        z <- data$cor
        xtitle <- "delta PSI"
        ytitle <- "log2FC expr"
        legendtitle <- "Expr-PSI change"
        
        ymin <- ceiling(max(abs(y[y<0]), na.rm=TRUE)) * -1
        ymax <- ceiling(max(y[y>0], na.rm=TRUE))
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, fill=z), color="black", pch=21, size=2, alpha=1.0, stroke=0.1) +
            geom_vline(xintercept=c(psi.delta*-1, psi.delta), linetype="dashed", color="black", size=0.2) +
            geom_hline(yintercept=c(gene.log2fc*-1, gene.log2fc), linetype="dashed", color="black", size=0.2) +
            scale_y_continuous(breaks=seq(ymin, ymax), limits=c(ymin, ymax)) +
            scale_fill_manual(values=colors) +
            labs(title=NULL, x=xtitle, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                plot.title = element_text(size=12, hjust=0.5),
                axis.line=element_line(colour = "black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                
    } else if(anno==TRUE){
        
        # Create labels
        df$label <- ifelse(df$cor != "Iso-Switch", paste(df$gene_short_name, df$event_type, sep="."), "")
        
        # Definition
        data <- df
        x <- data$mean.diff
        y <- data$log2fc.gene
        z <- data$cor
        xtitle <- "delta PSI"
        ytitle <- "log2FC expr"
        legendtitle <- "Expr-PSI change"
        label <- df$label
        
        ymin <- ceiling(max(abs(y[y<0]), na.rm=TRUE)) * -1
        ymax <- ceiling(max(y[y>0], na.rm=TRUE))
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, fill=z), color="black", pch=21, size=2, alpha=1.0, stroke=0.1) +
            ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 1.0, size=2, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
            geom_vline(xintercept=c(psi.delta*-1, psi.delta), linetype="dashed", color="black", size=0.2) +
            geom_hline(yintercept=c(gene.log2fc*-1, gene.log2fc), linetype="dashed", color="black", size=0.2) +
            scale_y_continuous(breaks=seq(ymin, ymax), limits=c(ymin, ymax)) +
            scale_fill_manual(values=colors) +
            labs(title=NULL, x=xtitle, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                plot.title = element_text(size=12, hjust=0.5),
                axis.line=element_line(colour = "black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                    
    }
    
    ##########################################
    
    # Save to new slots
    MarvelObject$DE$Cor$PSIvsExpr$Plot <- plot
    
    return(MarvelObject)
        
}
