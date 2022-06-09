#' @title Annotate reduced dimension space with gene expression values
#'
#' @description Annotates reduced dimension space, e.g., UMAP and tSNE, with gene expression values. Values will be automatically be log2-transformed prior to plotting.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param gene_short_name Character string. Gene name whose expression will be plotting.
#' @param point.size Numeric value. Size of data points. Default is \code{1}.
#' @param color.gradient Vector of character strings. Colors to indicate low, moderate, and high expression. Default is \code{c("grey90","blue","red")}.
#' @param type Character string. Type of reduced dimension space. Options are \code{"umap"} and \code{"tsne"}.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$adhocPlot$PCA$Gene}.
#'
#' @importFrom plyr join
#' @import ggplot2
#'
#' @export

PlotValues.PCA.Gene.10x <- function(MarvelObject, gene_short_name, point.size=0.1, color.gradient=c("grey90","blue","red"), type) {

    # Example aruguments
    MarvelObject <- MarvelObject
    df.coord <- MarvelObject$pca
    df.gene.norm <- MarvelObject$gene.norm.matrix
    gene_short_name <- gene_short_name
    point.size <- point.size
    color.gradient <- color.gradient
    type <- type
    
    # Example aruguments
    #MarvelObject <- marvel
    #df.coord <- MarvelObject$pca
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #gene_short_name <- "TPM2"
    #point.size <- 0.1
    #color.gradient <- c("grey","cyan","green","yellow","red")
    #type <- "tsne"
    
    ##########################################################################
    ###################### ANNOTATE EXPRESSION ###############################
    ##########################################################################
    
    # Subset relevant SJ
    df.gene.norm <- df.gene.norm[gene_short_name, ]
    
    # Tabulate counts
    df.gene.norm  <- data.frame("cell.id"=names(df.gene.norm),
                                "expr.gene.norm"=as.numeric(df.gene.norm),
                                stringsAsFactors=FALSE
                                )
    
    # Annotate coordinate table
    df.coord <- join(df.coord, df.gene.norm, by="cell.id", type="left")
    

    ##########################################################################
    ####################### SET THRESHOLD FOR COUNTS #########################
    ##########################################################################
    
    # Define margins before removing cells
    xmin <- min(df.coord$x); xmax <- max(df.coord$x)
    ymin <- min(df.coord$y); ymax <- max(df.coord$y)
    
    # Remove unexpressed genes
    #coord.xy.small <- coord.xy.small[which(coord.xy.small$count.exp >= 1), ]

    #if(min.gene.count >= 2) {
    
        # Remove lowly-expressed genes
        #coord.xy.small <- coord.xy.small[which(coord.xy.small$count.exp >= min.gene.count), ]
        
    #}
    
    # Censor lowly expressed gene
    #df.coord$expr.gene.norm[which(df.coord$gene.count < min.gene.count)] <- NA
    
    # Remove non-expressing/lowly expressed gene
    #df.coord <- df.coord[!is.na(df.coord$expr.gene.norm), ]
        
    ##########################################################################
    ################################ PLOT ####################################
    ##########################################################################
   
    # Reorder by expression
    df.coord <- df.coord[order(df.coord$expr.gene.norm), ]
        
    # Definitions
    data <- df.coord
    x <- data$x
    y <- data$y
    z <- log2(data$expr.gene.norm + 1)
    
    if(type=="umap"){
        
        xtitle <- "UMAP-1"
        ytitle <- "UMAP-2"
        
    } else if(type=="tsne"){
    
        xtitle <- "tSNE-1"
        ytitle <- "tSNE-2"
    
    }
    
    legendtitle <- "log2(expr)"
    
    # Plot
    plot <- ggplot() +
        geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size) +
        scale_color_gradientn(colors=color.gradient) +
        scale_x_continuous(limits=c(xmin, xmax)) +
        scale_y_continuous(limits=c(ymin, ymax)) +
        labs(title=NULL, x=xtitle, y=ytitle, color=legendtitle) +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            plot.title=element_text(hjust = 0.5, size=15),
            plot.subtitle=element_text(hjust = 0.5, size=15),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_text(size=12),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            #legend.position="none",
            legend.title=element_text(size=8),
            legend.text=element_text(size=8)
            )

    ##########################################################################

    # Save into new slot
    MarvelObject$adhocPlot$PCA$Gene <- plot

    # Return final object
    return(MarvelObject)

}
