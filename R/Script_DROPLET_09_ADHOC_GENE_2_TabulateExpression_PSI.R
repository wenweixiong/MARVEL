#' @title Dotplot of splice junction expression values for a specified gene
#'
#' @description Creates a dotplot of splice junction expression value of a specified gene across different cell groups. The gene and cell groups were defined earlier in \code{adhocGene.TabulateExpression.Gene.10x} function.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{adhocGene.TabulateExpression.Gene.10x} function.
#' @param min.pct.cells Numeric value. Percentage of cell expressing the splice junction in a cell group, below which, the value be re-coded as missing and appear will be omitted from the plot. A splice junction is considered to be expressed in a given cell if it has count >=1.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocGene$Expression$PSI$Table} and  \code{MarvelObject$adhocGene$Expression$PSI$Plot}.
#'
#' @importFrom plyr join
#' @importFrom stats aggregate
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
#' # SJ usage profiling
#' marvel.demo.10x <- adhocGene.TabulateExpression.PSI.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         min.pct.cells=10
#'                         )
#'
#' # Check output
#' marvel.demo.10x$adhocGene$Expression$PSI$Plot

adhocGene.TabulateExpression.PSI.10x <- function(MarvelObject, min.pct.cells=10) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.group.list <- MarvelObject$adhocGene$cell.group.list
    gene_short_name <- MarvelObject$adhocGene$gene_short_name
    min.pct.cells <- min.pct.cells
    
    # Example arguments
    #MarvelObject <- marvel
    #sj.metadata <- MarvelObject$sj.metadata
    #df.gene.count <- MarvelObject$gene.count.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #cell.group.list <- MarvelObject$adhocGene$cell.group.list
    #gene_short_name <- MarvelObject$adhocGene$gene_short_name
    #min.pct.cells <- 10
    
    ##########################################################################
    ########################### TABULATE SJ COUNTS ###########################
    ##########################################################################
    
    # Define SJ
    coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start==gene_short_name), "coord.intron"]
    
    # Subset relevant SJ
    df.sj.count <- df.sj.count[coord.introns, ]
    
    # Compute counts, % expr for each cell group
    .list <- list()
    
    for(i in 1:length(cell.group.list)) {
        
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.list[[i]]]
        
        # Compute counts
        sj.count.total <- apply(df.sj.count.small, 1, function(x) {sum(x)})
        
        # % cells expressing SJ
        n.cells.total <- ncol(df.sj.count.small)
        n.cells.expr.sj <- apply(df.sj.count.small, 1, function(x) {sum(x != 0)})
        pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
        
        # Save into data frame
        results <- data.frame("group"=names(cell.group.list)[i],
                              "coord.intron"=row.names(df.sj.count.small),
                              "n.cells.total"=n.cells.total,
                              "n.cells.expr.sj"=n.cells.expr.sj,
                              "pct.cells.expr.sj"=pct.cells.expr.sj,
                              "sj.count.total"=sj.count.total,
                              stringsAsFactors=FALSE
                              )
                              
        row.names(results) <- NULL
        
        # Save into list
        .list[[i]] <- results
        
    }
    
    results.sj.count <- do.call(rbind.data.frame, .list)
    
    ##########################################################################
    ########################## TABULATE GENE COUNTS ##########################
    ##########################################################################
        
    # Subset relevant gene
    df.gene.count <- df.gene.count[gene_short_name, , drop=FALSE]
    
    # Compute counts, % expr for each cell group
    .list <- list()
    
    for(i in 1:length(cell.group.list)) {
        
        # Subset cells
        df.gene.count.small <- df.gene.count[, cell.group.list[[i]], drop=FALSE]
        
        # Compute counts
        gene.count.total <- apply(df.gene.count.small, 1, function(x) {sum(x)})
        
        # % cells expressing SJ
        n.cells.total <- ncol(df.gene.count.small)
        n.cells.expr.gene <- apply(df.gene.count.small, 1, function(x) {sum(x != 0)})
        pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
        
        # Save into data frame
        results <- data.frame("group"=names(cell.group.list)[i],
                              "gene_short_name"=row.names(df.gene.count.small),
                              "n.cells.total"=n.cells.total,
                              "n.cells.expr.gene"=n.cells.expr.gene,
                              "pct.cells.expr.gene"=pct.cells.expr.gene,
                              "gene.count.total"=gene.count.total,
                              stringsAsFactors=FALSE
                              )
                              
        row.names(results) <- NULL
        
        # Save into list
        .list[[i]] <- results
        
    }
    
    results.gene.count <- do.call(rbind.data.frame, .list)
    
    ##########################################################################
    ################################# MERGE ##################################
    ##########################################################################
    
    # Annotate gene count
    results <- join(results.sj.count, results.gene.count[,c("group", "gene.count.total")], by="group", type="left")
    
    # Compute PSI
    results$psi <- round(results$sj.count.total/results$gene.count.total * 100, digits=2)
    
    # Remove lowly expressed SJ
    results$psi[which(results$pct.cells.expr.sj < min.pct.cells)] <- NA
    results <- results[!is.na(results$psi), ]
    
    # Re-order by % expressing cells
        # Reorder by % expressing cells
        #. <- aggregate(results$pct.cells.expr.sj, list(results$coord.intron), function(x) {sum(!is.na(x))})
        . <- aggregate(results$n.cells.expr.sj, list(results$coord.intron), function(x) {sum(x)})
        . <- .[order(.[,2], decreasing=TRUE), ]
        results$coord.intron <- factor(results$coord.intron, levels=.[,1])
                        
        # Reorder by cell group
        results$group <- factor(results$group, levels=names(cell.group.list))
        results <- results[order(results$group, results$coord.intron), ]
        
        # Annotate column no.
        coord.introns <- levels(results$coord.intron)
        . <- data.frame("coord.intron"=coord.introns,
                        "figure.column"=paste("SJ-", c(1:length(coord.introns)), sep=""),
                        stringsAsFactors=FALSE
                        )
        results <- join(results, ., by="coord.intron", type="left")
        
        # Reorder columns
        cols.1 <- c("group", "figure.column", "coord.intron")
        cols.2 <- setdiff(names(results), cols.1)
        results <- results[, c(cols.1, cols.2)]
            
    # Bubble plot
        # Definition
        data <- results
        x <- data$coord.intron
        y <- data$group
        z1 <- data$pct.cells.expr.sj
        z2 <- data$psi
        maintitle <- gene_short_name
        xtitle <- ""
        ytitle <- ""
        legendtitle.size <- "% cells"
        legendtitle.color <- "PSI"
        labels.y <- data$group
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, size=z1, color=z2)) +
            scale_color_gradientn(colors=c("gray","cyan","green","yellow","red")) +
            scale_y_discrete(limits = rev(levels(y))) +
            labs(title=maintitle, x=xtitle, y=ytitle, size=legendtitle.size, color=legendtitle.color) +
            theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                plot.title = element_text(size=12, hjust=0.5),
                axis.line.x=element_blank(),
                axis.line.y=element_line(colour = "black"),
                axis.ticks.x=element_blank(),
                axis.title=element_text(size=12),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=10),
                legend.text=element_text(size=8),
                legend.key=element_blank()
                ) +
            scale_size_area(breaks=c(25,50,75,100), limits=c(1, 100))
    
    ##################################################################
    
    # Save into new slots
    MarvelObject$adhocGene$Expression$PSI$Table <- results
    MarvelObject$adhocGene$Expression$PSI$Plot <- plot
    
    # Return final object
    return(MarvelObject)
            
}


