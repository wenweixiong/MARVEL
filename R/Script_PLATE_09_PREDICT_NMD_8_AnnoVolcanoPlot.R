#' @title Annotate volcano plot with nonsense-mediated decay (NMD) genes
#'
#' @description Annotate volcano plot generated from differential gene expression analysis with genes predicted to undergo splicing-induced NMD.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareExpr} function.
#' @param anno Logical value. If set to \code{TRUE}, selected gene names will be annotated on the plot as defined in \code{gene.label.x.below} and \code{gene.label.y.above}.
#' @param anno.gene_short_name Vector of character strings. When \code{anno} set to \code{TRUE}, the gene names to annotate on the plot.
#' @param label.size Numeric value. When \code{anno} set to \code{TRUE}, the size of gene labels.
#' @param point.size Numeric value. Size of data points. Default value is \code{1}.
#' @param xlabel.size Numeric value. Font size of the xtick labels. Default is \code{8}.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$NMD$AnnoVolcanoPlot$Table} and \code{MarvelObject$NMD$AnnoVolcanoPlot$Plot}.
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
#' marvel.demo <- AnnoVolcanoPlot(MarvelObject=marvel.demo)
#'
#' # Check outputs
#' head(marvel.demo$NMD$AnnoVolcanoPlot$Table)
#' marvel.demo$NMD$AnnoVolcanoPlot$Plot

AnnoVolcanoPlot <- function(MarvelObject, anno=FALSE, anno.gene_short_name=NULL, label.size=NULL, point.size=1.0, xlabel.size=8) {

    # Define arguments
    df <- MarvelObject$NMD$NMD.Expr$Table
    de <- MarvelObject$DE$Exp.Spliced$Table
    df.feature <- MarvelObject$GeneFeature
    anno <- anno
    anno.gene_short_name <- anno.gene_short_name
    point.size <- point.size
    xlabel.size <- xlabel.size
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$NMD$NMD.Expr$Table
    #de <- MarvelObject$DE$Exp.Spliced$Table
    #df.feature <- MarvelObject$GeneFeature
    #anno <- TRUE
    #anno.gene_short_name <- gene_short_names
    #label.size <- 2.5
    #xlabel.size <- 8
    #point.size <- 1
   
    # Subset NMD genes only
    df <- df[which(df$NMD==TRUE), ]
    
    # Convert factor to character for selected variables
    df$event_type <- as.character(df$event_type)
    df$NMD <- as.character(df$NMD)
    
    # Collapse genes with multiple events
        # Check for multi-events
        . <- as.data.frame(table(df$gene_id))
        gene_ids <- as.character(.[which(.$Freq >=2), "Var1"])
        
        # Collapse
        if(length(gene_ids) >= 1) {
            
            df$event_type[which(df$gene_id %in% gene_ids)] <- "Multi"
            df <- unique(df)
            
        }
    
    # Annotate DE gene result table
    de <- join(de, df[,c("gene_id", "NMD", "event_type")], by="gene_id", type="left")
    de$NMD[is.na(de$NMD)] <- "NOS"
    de$NMD[which(de$NMD==TRUE)] <- de$event_type[which(de$NMD==TRUE)]
    
    # Set factor levels
    levels <- intersect(c("SE", "RI", "A5SS", "A3SS", "Multi", "NOS"), unique(de$NMD))
    de$NMD <- factor(de$NMD, levels=levels)
    de <- de[order(de$NMD, decreasing=TRUE), ]
    
    # Set color scheme
    . <- data.frame("NMD"=levels, color=NA)
    .$color[which(.$NMD == "NOS")] <- "gray90"
    .$color[which(.$NMD != "NOS")] <- hue_pal()(length(levels)-1)
    colors <- .$color
    
    if(anno==FALSE) {
        
        # Boxplot
            # Definition
            data <- de
            x <- data$log2fc
            y <- -log10(data$p.val.adj)
            z <- data$NMD
            maintitle <- ""
            ytitle <- "-log10(p-value)"
            xtitle <- "log2FC"
            legendtitle <- "Event Type"
            
            xmin <- floor(min(x)) ; xmax <- ceiling(max(x)) ; xinterval <- 2
            
            if((xmin %% 2) != 0) {
                
                xmin <- xmin - 1 ; xmax <- xmax + 1
                
            }
            
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size) +
                geom_vline(xintercept=0, linetype="dashed", color="black", size=0.25) +
                scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                scale_color_manual(values=colors) +
                labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
                theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border=element_blank(),
                    plot.title=element_text(hjust = 0.5, size=12),
                    plot.subtitle=element_text(hjust = 0.5, size=12),
                    axis.line.y.left = element_line(color="black"),
                    axis.line.x = element_line(color="black"),
                    axis.title=element_text(size=12),
                    axis.text.x=element_text(size=xlabel.size, colour="black"),
                    axis.text.y=element_text(size=10, colour="black"),
                    legend.title=element_text(size=8),
                    legend.text=element_text(size=8)
                    )
                    
    } else {
        
        # Indicate labels
        de$label <- ifelse(de$gene_short_name %in% anno.gene_short_name, de$gene_short_name, "")
        
        # Boxplot
            # Definition
            data <- de
            x <- data$log2fc
            y <- -log10(data$p.val.adj)
            z <- data$NMD
            maintitle <- ""
            ytitle <- "-log10(p-value)"
            xtitle <- "log2FC"
            legendtitle <- "Event Type"
            label <- data$label
           
            xmin <- floor(min(x)) ; xmax <- ceiling(max(x)) ; xinterval <- 2
           
            if((xmin %% 2) != 0) {
               
               xmin <- xmin - 1 ; xmax <- xmax + 1
               
            }
           
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size) +
                geom_vline(xintercept=0, linetype="dashed", color="black", size=0.25) +
                ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0.01) +
                scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                scale_color_manual(values=colors) +
                labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
                theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border=element_blank(),
                    plot.title=element_text(hjust = 0.5, size=12),
                    plot.subtitle=element_text(hjust = 0.5, size=12),
                    axis.line.y.left = element_line(color="black"),
                    axis.line.x = element_line(color="black"),
                    axis.title=element_text(size=12),
                    axis.text.x=element_text(size=xlabel.size, colour="black"),
                    axis.text.y=element_text(size=10, colour="black"),
                    legend.title=element_text(size=8),
                    legend.text=element_text(size=8)
                    )
                    
    }
    
    # Subset NMD from DE results table
    de <- de[which(de$NMD != "NOS"), ]
    
    # Save to new slot
    MarvelObject$NMD$AnnoVolcanoPlot$Table <- de
    MarvelObject$NMD$AnnoVolcanoPlot$Plot <- plot
    
    # Return MARVEL object
    return(MarvelObject)
    
}
        
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
  
