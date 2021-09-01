#' @title Annotate Volcano Plot with Nonsense-Mediated Decay (NMD) Genes
#'
#' @description
#' \code{AnnoVolcanoPlot} plots annotates volcano plot generated from differential gene expression analysis with genes predicted to undergo splicing-induced NMD.
#'
#' @details
#' This function plots and annotates volcano plot generated from differential gene expression analysis with genes predicted to undergo splicing-induced NMD. Helpful to look for NMD genes based on degree of log2(fold change) and adjusted p-values. NMD genes with large log2(fold change) and small adjusted p-values may be of interest.
#'
#' @param MarvelObject S3 object generated from \code{CompareExpr} function.
#' @param anno Logical value. If set to \code{TRUE}, selected gene names will be annotated on the plot as defined in \code{gene.label.x.below} and \code{gene.label.y.above}.
#' @param gene.label.x.below Numeric value. Only applicable when \code{anno} set to \code{TRUE}. To indicate log2(fold change) value below which the NMD genes are annotated on the plot. Specified together with \code{gene.label.y.above}.
#' @param gene.label.y.above Numeric value. Only applicable when \code{anno} set to \code{TRUE}. To indicate -log10(p-value) above which the NMD genes are annotated on the plot. Specified together with \code{gene.label.x.below}.
#' @param gene.label.size Numeric value. Indicate size of gene labels specified using \code{gene.label.x.below} and \code{gene.label.y.above}.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$NMD$AnnoVolcanoPlot}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import ggplot2
#' @import ggrepel
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- AnnoVolcanoPlot(MarvelObject=marvel)
#'
#' # Check output
#' marvel$NMD$AnnoVolcanoPlot$Table[1:5, ]
#' marvel$NMD$AnnoVolcanoPlot$Plot

AnnoVolcanoPlot <- function(MarvelObject, anno=FALSE, gene.label.x.below=NULL, gene.label.y.above=NULL, gene.label.size=NULL) {

    # Define arguments
    df <- MarvelObject$NMD$NMD.Expr$Table
    de <- MarvelObject$DE$Exp$Table
    df.feature <- MarvelObject$GeneFeature
    anno <- anno
    gene.label.x.below <- gene.label.x.below
    gene.label.y.above <- gene.label.y.above
    
    # Example arguments
    #df <- marvel$NMD$NMD.Expr$Table
    #de <- marvel$DE$Exp$Table
    #df.feature <- marvel$GeneFeature
    #anno <- TRUE
    #gene.label.x.below <- 0
    #gene.label.y.above <- 10
    #gene.label.size <- 2.5
    
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
            
            xmin <- floor(min(x)) ; xmax <- ceiling(max(x)) ; xinterval <- 1
            
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z)) +
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
                    axis.text.x=element_text(size=10, colour="black"),
                    axis.text.y=element_text(size=10, colour="black"),
                    legend.title=element_text(size=8),
                    legend.text=element_text(size=8)
                    )
                    
    } else {
        
        # Indicate labels
        index <- which(de$log2fc < gene.label.x.below & -log10(de$p.val.adj) > gene.label.y.above & de$NMD != "NOS")
        de$label <- NA
        de$label[index] <- de$gene_short_name[index]
        
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
           
            xmin <- floor(min(x)) ; xmax <- ceiling(max(x)) ; xinterval <- 1
           
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z)) +
                geom_vline(xintercept=0, linetype="dashed", color="black", size=0.25) +
                geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=gene.label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0.01) +
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
                    axis.text.x=element_text(size=10, colour="black"),
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
        
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
  
