#' @title Heatmap of pseudo-bulk PSI values across different cell groups
#'
#' @description Heatmap of pseudo-bulk PSI values. x-axis represent main cell group, e.g. donor. y-axis represent sub-cell group, e.g. cell type
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param coord.intron Character string. Coordinates of splice junction whose expression will be plotted.
#' @param x.column Character string. Column name in \code{MarvelObject$sample.metadata} to indicate main-cell groups.
#' @param y.column Character string. Column name in \code{MarvelObject$sample.metadata} to indicate sub-cell groups.
#' @param x.levels Vector of character strings. Main cell groups to include for plotting.
#' @param y.levels Vector of character strings. Sub-cell groups to include for plotting.
#' @param xlabels.size Numeric value. Font size of heatmap column (main cell group) labels. Default is \code{8}.
#' @param ylabels.size Numeric value. Font size of heatmap row (sub-cell group) labels. Default is \code{8}.
#' @param min.pct.cells.gene.expr Numeric value. Percentage of cell expressing the gene in a pseudobulk, below which, the value be re-coded as missing and appear as grey on the heatmap.
#' @param min.n.cells Numeric value. Minimum number of cells in a cell group, below which, the cell group will be censored and appear as grey on the heatmap.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocPlot$Heatmap$Pseudobulk$Plot} and \code{MarvelObject$adhocPlot$Heatmap$Pseudobulk$Data}.
#'
#' @importFrom plyr join
#' @importFrom Matrix colSums
#' @import ggplot2
#' @import reshape2
#' @import pheatmap
#'
#' @export

PlotValues.PSI.Pseudobulk.Heatmap.10x <- function(MarvelObject, coord.intron, x.column, y.column, x.levels, y.levels, xlabels.size=8, ylabels.size=8, min.pct.cells.gene.expr=10, min.n.cells=10) {

    # Example arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    sj.metadata <- MarvelObject$sj.metadata
    df.sj.count <- MarvelObject$sj.count.matrix
    df.gene.count <- MarvelObject$gene.count.matrix
    coord.intron <- coord.intron
    y.column <- y.column
    x.column <- x.column
    y.levels <- y.levels
    x.levels <- x.levels
    ylabels.size <- ylabels.size
    xlabels.size <- xlabels.size
    min.pct.cells.gene.expr <- min.pct.cells.gene.expr
    min.n.cells <- min.n.cells
    
    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #sj.metadata <- MarvelObject$sj.metadata
    #df.sj.count <- MarvelObject$sj.count.matrix
    #df.gene.count <- MarvelObject$gene.count.matrix
    #coord.intron <- coord.intron
    #y.column <- y.column
    #x.column <- x.column
    #y.levels <- y.levels
    #x.levels <- x.levels
    #ylabels.size <- 8
    #xlabels.size <- 8
    #min.pct.cells.gene.expr <- 10
    
    ##########################################################################
    
    # Compute SJ counts
    df.sj.count <- df.sj.count[coord.intron, , drop=FALSE]
    sj.counts <- colSums(df.sj.count)
    
    # Compute gene counts
    gene_short_name <- unique(sj.metadata[which(sj.metadata$coord.intron %in% coord.intron), "gene_short_name.start"])
    
    df.gene.count <- df.gene.count[gene_short_name, , drop=FALSE]
    gene.counts <- colSums(df.gene.count)
    
    # Compute PSI for each pseudo-bulk
    .list.psi. <- list()
    
    for(j in 1:length(x.levels)) {
        
        # Define cell x.level
        x.level <- x.levels[j]
        
        .list.psi <- list()
        
        for(i in 1:length(y.levels)) {
            
            # Define cell y.level
            y.level <- y.levels[i]
            
            # Retrieve cell ids
            index.x <- which(sample.metadata[[x.column]] == x.level)
            index.y <- which(sample.metadata[[y.column]] == y.level)
            index <- intersect(index.x, index.y)
            cell.ids <- sample.metadata[index, "cell.id"]
            
            # Compute total cells
            n.cells <- length(cell.ids)
                            
            # Compute total SJ counts
            sj.counts.small <- sj.counts[which(names(sj.counts) %in% cell.ids)]
            sj.counts.total <- sum(sj.counts.small)
            n.cells.sj.expr <- sum(sj.counts.small != 0)
            pct.cells.sj.expr <- round(n.cells.sj.expr/n.cells * 100, digits=1)
                
            # Compute total gene counts
            gene.counts.small <- gene.counts[which(names(gene.counts) %in% cell.ids)]
            gene.counts.total <- sum(gene.counts.small)
            n.cells.gene.expr <- sum(gene.counts.small != 0)
            pct.cells.gene.expr <- round(n.cells.gene.expr/n.cells * 100, digits=1)
                
            # Compute PSI
            psi <- sj.counts.total/gene.counts.total * 100
            
            # Censor PSI if gene not expressed
            index.l <- pct.cells.gene.expr < min.pct.cells.gene.expr | is.na(pct.cells.gene.expr) | n.cells < min.n.cells
            
            if(index.l) {
                
                psi <- NA
                
            }
                
            # Tabulate results
            results <- data.frame("x.level"=x.level,
                                  "y.level"=y.level,
                                  "psi"=psi,
                                  stringsAsFactors=FALSE
                                  )
                
            # Save results
            .list.psi[[i]] <- results
            
        }
                    
        # Collapse into data frame
        results. <- do.call(rbind.data.frame, .list.psi)
            
        # Save into list
        .list.psi.[[j]] <- results.
            
    }
    
    # Merge
    results <- do.call(rbind.data.frame, .list.psi.)
    
    # Set factor levels
    results$x.level <- factor(results$x.level, levels=x.levels)
    results$y.level <- factor(results$y.level, levels=y.levels)

    # Convert to heatmap matrix
    results <- dcast(data=results, formula=y.level ~ x.level, value.var="psi")
    row.names(results) <- results$y.level
    results$y.level <- NULL
    #results <- as.matrix(results)
    
    # Heatmap
    plot <- pheatmap(results,
                     cluster_rows=FALSE,
                     cluster_cols=FALSE,
                     scale="row",
                     fontsize_row=ylabels.size,
                     fontsize_col=xlabels.size,
                     border_color="white",
                     legend=TRUE,
                     show_rownames=TRUE,
                     angle_col=90,
                     silent=TRUE
                     #color=myColor,
                     #breaks=myBreaks
                     )
    
    ##########################################################################

    # Save into new slots
    MarvelObject$adhocPlot$Heatmap$Pseudobulk$Plot <- plot
    MarvelObject$adhocPlot$Heatmap$Pseudobulk$Data <- results
    
    # Return final object
    return(MarvelObject)

}
