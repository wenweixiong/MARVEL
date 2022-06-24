#' @title Violin plot of single-cell PSI values vs cell groups
#'
#' @description Violin plot of single-cell PSI values on the y-axis against cell groups on the x-axis.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param cell.group.list List of character strings. Each element in the list is a cell group. The label of each element of the list represents the cell group label.
#' @param gene_short_name Character string. Gene name whose expression will be plotted.
#' @param cell.group.colors Vector of character strings. Colors of cell groups and should be same length as \code{cell.group.list}. Default \code{ggplot2} colors are used.
#' @param xlabels.size Numeric value. Font size of x-tick labels. Default is \code{10}.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocPlot$ViolinPlot$SingleCell$Plot} and \code{MarvelObject$adhocPlot$ViolinPlot$SingleCell$Data}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @importFrom grDevices hcl
#'
#' @export

PlotValues.Gene.SingleCell.10x <- function(MarvelObject, cell.group.list, gene_short_name, cell.group.colors=NULL, xlabels.size=8, alpha=0.5) {

    # Example arguments
    MarvelObject <- MarvelObject
    df.gene.norm <- MarvelObject$gene.norm.matrix
    cell.group.list <- cell.group.list
    gene_short_name <- gene_short_name
    cell.group.colors <- cell.group.colors
    xlabels.size <- xlabels.size
    alpha <- alpha
    
    # Example arguments
    #MarvelObject <- marvel
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #cell.group.list <- group.list
    #gene_short_name <- "EZH2"
    #cell.group.colors <- NULL
    #xlabels.size <- 8
    
    ##########################################################################
    
    # Retrieve cell groups
    cell.groups <- names(cell.group.list)
    
    .list <- list()
    
    for(i in 1:length(cell.groups)) {
        
        # Define cell group
        cell.group <- cell.groups[i]
        
        # Retrieve cell ids
        cell.ids <- cell.group.list[[cell.group]]
        
        # Retrieve gene expression
        . <- df.gene.norm[gene_short_name, cell.ids]
        results <- data.frame("cell.group"=cell.group,
                              "cell.id"=names(.),
                              "exp"=as.numeric(.),
                              stringsAsFactors=FALSE
                              )
        
        # Save into list
        .list[[i]] <- results

    }
    
    results <- do.call(rbind.data.frame, .list)
    
    # Set factor levels
    results$cell.group <- factor(results$cell.group, levels=names(cell.group.list))
    
    # Create x-labels
        # Tabulate total cells
        n.cell.total <- as.data.frame(table(results$cell.group))[["Freq"]]
        
        # Tabulate expressing cells
        . <- by(results$exp, results$cell.group, function(x) {sum(x != 0)})
        n.cells.expr <- as.numeric(.)
        pct.cell.expr <- round(n.cells.expr/n.cell.total * 100, digits=0)
        
        # Create labels
        xlabels <- paste(levels(results$cell.group), "\n", n.cells.expr, "/", n.cell.total, "(", pct.cell.expr, "%)", sep="")
        
    # Violin plot
        # Definition
        data <- results
        x <- data$cell.group
        y <- log2(data$exp + 1)
        z <- data$cell.group
        maintitle <- ""
        ytitle <- "log2(Norm. UMI counts + 1)"
        xtitle <- ""
        xlabels <- xlabels

        # Color scheme
        if(is.null(cell.group.colors[1])) {
        
            gg_color_hue <- function(n) {
              hues = seq(15, 375, length = n + 1)
              hcl(h = hues, l = 65, c = 100)[1:n]
            }
            
            n = length(levels(z))
            cols = gg_color_hue(n)
        
        } else {
            
            cols <- cell.group.colors
            
        }

        # Plot
        plot <- ggplot() +
            geom_violin(data, mapping=aes(x=x, y=y, fill=z), scale="width") +
            geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.3, height=0), size=0.1, alpha=alpha) +
            stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
            scale_fill_manual(values=cols) +
            scale_x_discrete(labels=xlabels) +
            #scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=12),
                plot.subtitle=element_text(hjust = 0.5, size=12),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=xlabels.size, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.position="none"
                )
    
    ##########################################################################

    # Save into new slots
    MarvelObject$adhocPlot$ViolinPlot$SingleCell$Plot <- plot
    MarvelObject$adhocPlot$ViolinPlot$SingleCell$Data <- results
    
    # Return final object
    return(MarvelObject)

}
