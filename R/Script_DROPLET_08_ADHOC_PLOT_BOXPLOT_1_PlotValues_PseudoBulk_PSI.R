#' @title Boxplot of pseudo-bulk PSI values vs cell groups
#'
#' @description Boxplot of pseudo-bulk PSI values on the y-axis against cell groups on the x-axis.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param cell.group.list List of character strings. Each element in the list (main list) is a list (sub-list). Each sub-list in turn is a list of pseudobulk samples represented by cell IDs. Each element in the main list is a cell group and the label of each element of the main list represents the cell group label.
#' @param coord.intron Character string. Coordinates of splice junction whose expression will be plotted.
#' @param cell.group.colors Vector of character strings. Colors of cell groups and should be same length as \code{cell.group.list}. Default \code{ggplot2} colors are used.
#' @param xlabels.size Numeric value. Font size of x-tick labels. Default is \code{10}.
#' @param min.pct.cells.gene.expr Numeric value. Percentage of cell expressing the gene in a pseudobulk, below which, the pseudobulk will be omitted from plotting. Default is \code{10}.
#' @param min.n.cells.gene.expr Numeric value. Number of cell expressing the gene in a pseudobulk, below which, the pseudobulk will be omitted from plotting. Default is \code{3}.
#' @param min.gene.counts.total Numeric value. Totol gene counts in a pseudobulk, below which, the pseudobulk will be omitted from plotting. Default is \code{3}.
#' @param method Character string. Statistical test for all possible pair-wise comparisons. Options are \code{"t.test"} (default) or \code{"wilcox"}.
#' @param p.adjust.method Character string. Method for multiple testing adjustment as per \code{method} option of \code{p.adjust} function. Default is \code{"fdr"}.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Plot}, \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Stats}, and \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Data}.
#'
#' @importFrom plyr join
#' @importFrom stats p.adjust p.adjust.methods pairwise.t.test pairwise.wilcox.test
#' @import ggplot2
#' @importFrom grDevices hcl
#' @import Matrix
#'
#' @export

PlotValues.PSI.Pseudobulk.10x <- function(MarvelObject, cell.group.list, coord.intron, cell.group.colors=NULL, xlabels.size=10, min.pct.cells.gene.expr=10, min.n.cells.gene.expr=3, min.gene.counts.total=3, method="t.test", p.adjust.method="fdr") {

    # Define arguments
    MarvelObject <- MarvelObject
    sj.metadata <- MarvelObject$sj.metadata
    df.sj.count <- MarvelObject$sj.count.matrix
    df.gene.count <- MarvelObject$gene.count.matrix
    cell.group.list <- cell.group.list
    coord.intron <- coord.intron
    cell.group.colors <- cell.group.colors
    xlabels.size <- xlabels.size
    min.pct.cells.gene.expr <- min.pct.cells.gene.expr
    min.n.cells.gene.expr <- min.n.cells.gene.expr
    min.gene.counts.total <- min.gene.counts.total
    method <- method
    p.adjust.method <- p.adjust.method
    
    # Example arguments
    #MarvelObject <- marvel
    #sj.metadata <- MarvelObject$sj.metadata
    #df.sj.count <- MarvelObject$sj.count.matrix
    #df.gene.count <- MarvelObject$gene.count.matrix
    #cell.group.list <- group.list
    #coord.intron <- "chr2:190975888:190976839"
    #cell.group.colors <- NULL
    #xlabels.size <- 10
    #min.pct.cells.gene.expr <- 10
    #min.n.cells.gene.expr <- 10
    #min.gene.counts.total <- 3
    #method <- "t.test"
    #p.adjust.method <- "none"
    
    ##########################################################################
    
    # Compute SJ counts
    df.sj.count <- df.sj.count[coord.intron, , drop=FALSE]
    sj.counts <- Matrix::colSums(df.sj.count)
    
    # Compute gene counts
    gene_short_name <- unique(sj.metadata[which(sj.metadata$coord.intron %in% coord.intron), "gene_short_name.start"])
    
    df.gene.count <- df.gene.count[gene_short_name, , drop=FALSE]
    gene.counts <- Matrix::colSums(df.gene.count)
    
    # Compute PSI for each pseudo-bulk
    .list.psi. <- list()
    
    for(j in 1:length(cell.group.list)) {
        
        # Retrieve cell group
        cell.group.list.small <- cell.group.list[[j]]
        
            .list.psi <- list()
            
            for(i in 1:length(cell.group.list.small)) {
                
                # Retrieve cell ids
                cell.ids <- cell.group.list.small[[i]]
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
                
                # Tabulate results
                results <- data.frame("sample.id"=names(cell.group.list.small)[i],
                                      "n.cells.total"=n.cells,
                                      "sj.counts.total"=sj.counts.total,
                                      "n.cells.sj.expr"=n.cells.sj.expr,
                                      "pct.cells.sj.expr"=pct.cells.sj.expr,
                                      "gene.counts.total"=gene.counts.total,
                                      "n.cells.gene.expr"=n.cells.gene.expr,
                                      "pct.cells.gene.expr"=pct.cells.gene.expr,
                                      "psi"=psi,
                                      stringsAsFactors=FALSE
                                      )
                
                # Save results
                .list.psi[[i]] <- results
                
            }
            
            # Collapse into data frame
            results. <- do.call(rbind.data.frame, .list.psi)
            
            # Indicate cell group
            . <- data.frame("cell.group"=names(cell.group.list[j]))
            results. <- cbind.data.frame(., results.)
            
            # Save into list
            .list.psi.[[j]] <- results.
            
        }
    
    # Merge
    results <- do.call(rbind.data.frame, .list.psi.)
    
    # Set factor levels
    levels <- names(cell.group.list)
    results$cell.group <- factor(results$cell.group, levels)
    
    # Censor cell groups not expressing gene
    index <- which(results$pct.cells.gene.expr < min.pct.cells.gene.expr | results$n.cells.gene.expr < min.n.cells.gene.expr | results$gene.counts.total < min.gene.counts.total)
    results$psi[index] <- NA
    
    # Compute sample size
    results.small <- results[!is.na(results$psi), ]
    . <- as.data.frame(table(results.small$cell.group))
    xlabels <- paste(.[,1], "\n(n=", .[,2], ")" , sep="")
    
    # Ensure there are at least two groups with at least one expressing cell
    n.cell.groups.expr <- sum(.[,2] != 0)
    
    if(n.cell.groups.expr >= 2) {

        # Boxplot
            # Definition
            data <- results
            x <- data$cell.group
            y <- data$psi
            z <- data$cell.group
            maintitle <- ""
            ytitle <- "Pseudobulk PSI (%)"
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
                geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), outlier.shape=NA) +
                geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=1) +
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
        
        # Statistical test
        if(method=="t.test"){
            
            stats <- pairwise.t.test(x=y, g=x, p.adjust.method=p.adjust.method)
            
        } else if(method=="wilcox"){
           
           stats <- pairwise.wilcox.test(x=y, g=x, p.adjust.method=p.adjust.method)
            
        }
        
        stats <- stats$p.value
        . <- data.frame("V1"=row.names(stats))
        stats <- cbind.data.frame(., stats)
        names(stats)[1] <- ""
        row.names(stats) <- NULL
        
        ##########################################################################

        # Save into new slots
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Plot <- plot
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Stats <- stats
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Data <- results
        
    } else {
        
        message("No cells expressing gene for plotting")
        
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Plot <- NULL
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Stats <- NULL
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$PSI$Data <- results
        
    }
    
    # Return final object
    return(MarvelObject)

}
