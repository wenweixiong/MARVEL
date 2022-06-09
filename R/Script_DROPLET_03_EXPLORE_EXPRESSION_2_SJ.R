#' @title Plot splice junction expression distribution
#'
#' @description Generates a plot of splice junction expression distribution (percentage of cells expressing a particular splice junction) to determine splice junction expression threshold for downstream differential splice junction analysis.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group) of downstream differential splice junction analysis.
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2 of downstream differential splice junction analysis.
#' @param min.pct.cells.genes Numeric value. Minimum percentage of cells in which the gene is expressed for that gene to be included for splice junction expression distribution analysis. Expressed genes defined as genes with non-zero normalised UMI counts. This threshold may be determined from \code{PlotPctExprCells.SJ.10x} function.
#' @param min.pct.cells.sj Numeric value. Minimum percentage of cells in which the splice junction is expressed for that splice junction to be included for splice junction expression distribution analysis. Expressed splice junctions defined as splice junctions with raw UMI counts >= 1.
#'
#' @return An object of class S3 with a new slots \code{MarvelObject$pct.cells.expr$SJ$Plot} and \code{MarvelObject$pct.cells.expr$SJ$Data}
#'
#' @importFrom plyr join
#' @import ggplot2
#'
#' @export

PlotPctExprCells.SJ.10x <- function(MarvelObject, cell.group.g1, cell.group.g2, min.pct.cells.genes=10, min.pct.cells.sj=10) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.group.g1 <- cell.group.g1
    cell.group.g2 <- cell.group.g2
    min.pct.cells.genes <- min.pct.cells.genes
    min.pct.cells.sj <- min.pct.cells.sj
    
    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #sj.metadata <- MarvelObject$sj.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #cell.group.g1 <- cell.ids.1
    #cell.group.g2 <- cell.ids.2
    #min.pct.cells.genes <- 10
    #min.pct.cells.sj <- 5
    
    ################################################################
    
    # Compute num. of cells in which gene is expressed: Group 1
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g1]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .
    
    # Compute num. of cells in which gene is expressed: Group 2
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g2]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g2),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
        
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Censor non-expressing genes
    results <- results[which(results$pct.cells.expr > min.pct.cells.genes), ]
    
    # Save as new object
    results.genes <- results
    
    ################################################################
    
    # Compute num. of cells in which SJ is expressed: Group 1
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.g1]
        
        # Subset expressed genes
        gene_short_names <- results.genes[which(results.genes$cell.group=="cell.group.g1"), "gene_short_name"]
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "coord.intron"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .

    # Compute num. of cells in which SJ is expressed: Group 2
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.g2]
        
        # Subset expressed genes
        gene_short_names <- results.genes[which(results.genes$cell.group=="cell.group.g2"), "gene_short_name"]
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "coord.intron"=names(.),
                        "n.cells.total"=length(cell.group.g2),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
            
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Censor non-expressing genes
    results <- results[which(results$pct.cells.expr > min.pct.cells.sj), ]
    
    # Save as new object
    results.sj <- results
        
    ################################################################
    
    # Density plot
        # Definitions
        data <- results.sj
        x <- data$pct.cells.expr
        z <- data$cell.group
        maintitle <- ""
        ytitle <- "Density"
        xtitle <- "SJ Expressed (% Cells)"
        xmin <- 0 ; xmax <- max(x) ; xinterval <- 10
        legendtitle <- ""
       
        # Plot
        plot <- ggplot() +
            geom_density(data, mapping=aes(x=x, color=z), alpha=0.5) +
            scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
            labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
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
                axis.text.x=element_text(size=8, colour="black"),
                axis.text.y=element_text(size=8, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                    
    # Save into new slot
    MarvelObject$pct.cells.expr$SJ$Plot <- plot
    MarvelObject$pct.cells.expr$SJ$Data <- results
    
    # Return final object
    return(MarvelObject)
            
}


