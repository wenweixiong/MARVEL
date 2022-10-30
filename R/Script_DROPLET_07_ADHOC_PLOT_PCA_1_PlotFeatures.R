#' @title Annotate reduced dimension space with cell feature
#'
#' @description Annotates reduced dimension space, e.g., UMAP and tSNE, with cell features such as cell group, donor ID, sample ID, etc.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param cell.group.list List of character strings. Each element of the list is a vector of cell IDs corresponding to a feature, e.g. cell group. The names of the element will be the cell feature label.
#' @param legendtitle Character string. Legend title. Default is \code{"Cell group"}.
#' @param alpha Numeric value. Transparency of the data points. Takes any values between 0-1 whereby 0 is totally transparent and 1 is opaque. Default is \code{0.75}.
#' @param point.size Numeric value. Size of data points. Default is \code{1}.
#' @param point.stroke Numeric value. Outline thickness of data points.  Default is \code{0.1}.
#' @param point.colors Vector of character strings. Colors of cell groups and should be same length as \code{cell.group.list}. Default \code{ggplot2} colors are used.
#' @param point.size.legend Numeric value. Size of legend keys. Default is \code{2}.
#' @param type Character string. Type of reduced dimension space. Options are \code{"umap"} and \code{"tsne"}.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$adhocPlot$PCA$CellGroup}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @importFrom grDevices hcl
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
#' # Define cell groups
#'     # Retrieve sample metadata
#'     sample.metadata <- marvel.demo.10x$sample.metadata
#'
#'     # iPSC
#'     index <- which(sample.metadata$cell.type=="iPSC")
#'     cell.ids.1 <- sample.metadata[index, "cell.id"]
#'     length(cell.ids.1)
#'
#'     # Cardio day 10
#'     index <- which(sample.metadata$cell.type=="Cardio day 10")
#'     cell.ids.2 <- sample.metadata[index, "cell.id"]
#'     length(cell.ids.2)
#'
#'     # Save into list
#'     cell.group.list <- list("iPSC"=cell.ids.1,
#'                             "Cardio d10"=cell.ids.2
#'                             )
#'
#' # Plot cell groups
#' marvel.demo.10x <- PlotValues.PCA.CellGroup.10x(
#'                             MarvelObject=marvel.demo.10x,
#'                             cell.group.list=cell.group.list,
#'                             legendtitle="Cell group",
#'                             type="tsne"
#'                             )
#'
#' # Check output
#' marvel.demo.10x$adhocPlot$PCA$CellGroup

PlotValues.PCA.CellGroup.10x <- function(MarvelObject, cell.group.list, legendtitle="Cell group", alpha=0.75, point.size=1.0, point.stroke=0.1, point.colors=NULL, point.size.legend=2, type) {

    # Example aruguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$pca
    cell.group.list <- cell.group.list
    legendtitle <- legendtitle
    alpha <- alpha
    point.size <- point.size
    point.stroke <- point.stroke
    point.colors <- point.colors
    point.size.legend <- point.size.legend
    type <- type
    
    # Example aruguments
    #MarvelObject <- marvel
    #df <- MarvelObject$pca
    #cell.group.list <- cell.group.list[1]
    #legendtitle <- "Cell type"
    #alpha <- 0.75
    #point.size <- 0.2
    #point.stroke <- 0
    #point.colors <- NULL
    #point.size.legend <- 2
    #type <- "tsne"
    
    ##########################################################################
    ########################### SUBSET CELLS #################################
    ##########################################################################
    
    # Create reference table for annotation
    .list <- list()
    
    for(i in 1:length(cell.group.list)) {
        
        . <- data.frame("cell.id"=cell.group.list[[i]],
                        "group"=names(cell.group.list)[i]
                        )
                        
        .list[[i]] <- .
        
    }
    
    ref <- do.call(rbind.data.frame, .list)
    
    # Annotate coordinate table
    df <- join(df, ref, by="cell.id", type="left")
        
    # Set factor levels
    df$group <- factor(df$group, levels=names(cell.group.list))
    
    # Report progress
    n.anno.missing <- sum(is.na(df$group))
    
    if(n.anno.missing==0) {
        
        message("All cells defined with coordinates found")
        
    } else {
        
        message(paste(n.anno.missing, " cells defined with no coordinates found", sep=""))
        
    }
    
    ##########################################################################

    # Definition
    data <- df
    x <- data$x
    y <- data$y
    z <- data$group
    
    if(type=="umap"){
        
        xtitle <- "UMAP-1"
        ytitle <- "UMAP-2"
        
    } else if(type=="tsne"){
    
        xtitle <- "tSNE-1"
        ytitle <- "tSNE-2"
    
    }
    
    # Define color scheme
    if(is.null(point.colors[1])) {
      
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
                            
      n = length(levels(z))
      point.colors = gg_color_hue(n)

    }

    # No. of columns for legend
    if(length(unique(z)) > 20) {

        ncol.legend <- 3

    } else if(length(unique(z)) > 10) {
        
        ncol.legend <- 2
    
    } else if(length(unique(z)) <= 10) {
        
        ncol.legend <- 1
        
    }
    
    # Plot (with legends)
    plot <- ggplot() +
          geom_point(data, mapping=aes(x=x, y=y, fill=z), color="black", pch=21, size=point.size, alpha=alpha, stroke=point.stroke) +
          scale_fill_manual(values=point.colors) +
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
              ) +
           guides(fill = guide_legend(override.aes=list(size=point.size.legend, alpha=alpha, stroke=point.stroke), ncol=ncol.legend))

    ##########################################################################

    # Save into new slot
    MarvelObject$adhocPlot$PCA$CellGroup <- plot

    # Return final object
    return(MarvelObject)

}
