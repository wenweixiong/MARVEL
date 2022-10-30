#' @title Boxplot of pseudo-bulk PSI values vs cell groups
#'
#' @description Boxplot of pseudo-bulk PSI values on the y-axis against cell groups on the x-axis.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param cell.group.list List of character strings. Each element in the list (main list) is a list (sub-list). Each sub-list in turn is a list of pseudobulk samples represented by cell IDs. Each element in the main list is a cell group and the label of each element of the main list represents the cell group label.
#' @param gene_short_name Character string. Name of gene to plot.
#' @param log2.transform Logical value. If set to \code{TRUE} (default), normalised gene expression values will be off-set by 1 and then log2-transformed prior to plotting.
#' @param cell.group.colors Vector of character strings. Colors of cell groups and should be same length as \code{cell.group.list}. Default \code{ggplot2} colors are used.
#' @param xlabels.size Numeric value. Font size of x-tick labels. Default is \code{10}.
#' @param min.n.cells.total Numeric value. Minimum number of cells in a pseudobulk, below which, the pseudobulk will be omitted from plotting. Default is \code{10}.
#' @param method Character string. Statistical test for all possible pair-wise comparisons. Options are \code{"t.test"} (default) or \code{"wilcox"}.
#' @param p.adjust.method Character string. Method for multiple testing adjustment as per \code{method} option of \code{p.adjust} function. Default is \code{"fdr"}.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Plot}, \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Stats}, and \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Data}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @importFrom grDevices hcl
#' @import Matrix
#'
#' @export

PlotValues.Gene.Pseudobulk.10x <- function(MarvelObject, cell.group.list, gene_short_name, log2.transform=TRUE, cell.group.colors=NULL, xlabels.size=10, min.n.cells.total=10, method="t.test", p.adjust.method="fdr") {

    # Define arguments
    MarvelObject <- MarvelObject
    df.gene.norm <- MarvelObject$gene.norm.matrix
    cell.group.list <- cell.group.list
    gene_short_name <- gene_short_name
    cell.group.colors <- NULL
    xlabels.size <- xlabels.size
    min.n.cells.total <-  min.n.cells.total
    p.adjust.method <- p.adjust.method
    log2.transform <- log2.transform
    
    # Example arguments
    #MarvelObject <- marvel
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #cell.group.list <- cell.group.list
    #gene_short_name <- "CEP135"
    #cell.group.colors <- NULL
    #xlabels.size <- 10
    #min.n.cells.total <- 25
    #p.adjust.method <- "none"
    
    ##########################################################################
    
    # Subset gene
    df.gene.norm <- df.gene.norm[gene_short_name, , drop=FALSE]
    
    # Transform values
    if(log2.transform==TRUE) {
        
        df.gene.norm <- log2(df.gene.norm + 1)
        
    }
    
    # Compute mean expression for each pseudo-bulk
    .list <- list()
    
    for(j in 1:length(cell.group.list)) {
        
        # Retrieve cell group
        cell.group.list.small <- cell.group.list[[j]]
        
            .list. <- list()
            
            for(i in 1:length(cell.group.list.small)) {
                
                # Retrieve cell ids
                cell.ids <- cell.group.list.small[[i]]
                n.cells <- length(cell.ids)
                
                # Compute average expression
                mean.expr.gene.norm <- mean(df.gene.norm[, cell.ids])
                
                # Tabulate results
                results <- data.frame("sample.id"=names(cell.group.list.small)[i],
                                      "n.cells.total"=n.cells,
                                      "mean.expr.gene.norm"=mean.expr.gene.norm,
                                      stringsAsFactors=FALSE
                                      )
                
                # Save results
                .list.[[i]] <- results
                
            }
            
            # Collapse into data frame
            results. <- do.call(rbind.data.frame, .list.)
            
            # Indicate cell group
            . <- data.frame("cell.group"=names(cell.group.list[j]))
            results. <- cbind.data.frame(., results.)
            
            # Save into list
            .list[[j]] <- results.
            
        }
    
    # Merge
    results <- do.call(rbind.data.frame, .list)
    
    # Set factor levels
    levels <- names(cell.group.list)
    results$cell.group <- factor(results$cell.group, levels)
    
    # Censor cell groups not expressing gene
    index <- which(results$n.cells.total < min.n.cells.total)
    results$mean.expr.gene.norm[index] <- NA
    
    # Compute sample size
    results.small <- results[!is.na(results$mean.expr.gene.norm), ]
    . <- as.data.frame(table(results.small$cell.group))
    xlabels <- paste(.[,1], "\n(n=", .[,2], ")" , sep="")
    
    # Ensure there are at least two groups with at least one expressing cell
    n.cell.groups.expr <- sum(.[,2] != 0)
    
    if(n.cell.groups.expr >= 2) {

        # Boxplot
            # Definition
            data <- results
            x <- data$cell.group
            y <- data$mean.expr.gene.norm
            z <- data$cell.group
            maintitle <- ""
            ytitle <- "Pseudobulk norm. expression"
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
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Plot <- plot
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Stats <- stats
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Data <- results
        
    } else {
        
        message("No cells expressing gene for plotting")
        
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Plot <- NULL
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Stats <- NULL
        MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Data <- results
        
    }
    
    # Return final object
    return(MarvelObject)

}
