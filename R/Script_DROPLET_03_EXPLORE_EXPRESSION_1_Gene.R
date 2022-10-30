#' @title Plot gene expression distribution
#'
#' @description Generates a plot of gene expression distribution (percentage of cells expressing a particular gene) to determine normalised gene expression threshold for downstream differential splice junction analysis.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group) of downstream differential splice junction analysis.
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2 of downstream differential splice junction analysis.
#' @param min.pct.cells Numeric value. Minimum percentage of cells in which the gene is expressed for that gene to be included for gene expression distribution analysis. Expressed genes defined as genes with non-zero normalised UMI counts.
#'
#' @return An object of class S3 with a new slots \code{MarvelObject$pct.cells.expr$Gene$Plot} and \code{MarvelObject$pct.cells.expr$Gene$Data}.
#'
#' @importFrom plyr join
#' @import ggplot2
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
#'     # Group 1 (reference)
#'     index <- which(sample.metadata$cell.type=="iPSC")
#'     cell.ids.1 <- sample.metadata[index, "cell.id"]
#'     length(cell.ids.1)
#'
#'     # Group 2
#'     index <- which(sample.metadata$cell.type=="Cardio day 10")
#'     cell.ids.2 <- sample.metadata[index, "cell.id"]
#'     length(cell.ids.2)
#'
#' # Explore % of cells expressing genes
#' marvel.demo.10x <- PlotPctExprCells.Genes.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         cell.group.g1=cell.ids.1,
#'                         cell.group.g2=cell.ids.2,
#'                         min.pct.cells=5
#'                         )
#'
#' # Check output
#' marvel.demo.10x $pct.cells.expr$Gene$Plot
#' head(marvel.demo.10x $pct.cells.expr$Gene$Data)

PlotPctExprCells.Genes.10x <- function(MarvelObject, cell.group.g1, cell.group.g2, min.pct.cells=1) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    cell.group.g1 <- cell.group.g1
    cell.group.g2 <- cell.group.g2
    min.pct.cells <- min.pct.cells
    
    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #cell.group.g1 <- cell.ids.1
    #cell.group.g2 <- cell.ids.2
    #min.pct.cells <- 10
    
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
    results <- results[which(results$pct.cells.expr > min.pct.cells), ]

    # Density plot
        # Definitions
        data <- results
        x <- data$pct.cells.expr
        z <- data$cell.group
        maintitle <- ""
        ytitle <- "Density"
        xtitle <- "Gene Expressed (% Cells)"
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
    MarvelObject$pct.cells.expr$Gene$Plot <- plot
    MarvelObject$pct.cells.expr$Gene$Data <- results
    
    # Return final object
    return(MarvelObject)
            
}


