#' @title Dotplot of gene expression values for a specified gene
#'
#' @description Creates a dotplot of average expression value of a specified gene across different cell groups.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param cell.group.list List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group.
#' @param gene_short_name Character string. Gene names whose expression will be plotted.
#' @param log2.transform Logical value. If set to \code{TRUE} (default), normalised gene expression values will be off-set by 1 and then log2-transformed prior to plotting.
#' @param min.pct.cells Numeric value. Percentage of cell expressing the gene in a cell group, below which, the value be re-coded as missing and appear will be omitted from the plot. A gene is considered to be expressed in a given cell if it has non-zero normalised count.
#' @param downsample Logical value. If set to \code{TRUE}, the number of cells in each cell group will be down-sampled so that all cell groups will have the same number of cells. The number of cells to down-sample will be based on the smallest cell group. Default is \code{FALSE}.
#' @param seed Numeric value. Random number generator to be fixed for down-sampling.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocGene$Expression$Gene$Table}, \code{MarvelObject$adhocGene$Expression$Gene$Plot}, \code{MarvelObject$adhocGene$cell.group.list}, and \code{MarvelObject$adhocGene$gene_short_name}.
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
#' # Gene expression profiling
#' marvel.demo.10x <- adhocGene.TabulateExpression.Gene.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         cell.group.list=cell.group.list,
#'                         gene_short_name="TPM2",
#'                         min.pct.cells=10,
#'                         downsample=TRUE
#'                         )
#'
#' # Check output
#' marvel.demo.10x$adhocGene$Expression$Gene$Plot
#' marvel.demo.10x$adhocGene$Expression$Gene$Table

adhocGene.TabulateExpression.Gene.10x <- function(MarvelObject, cell.group.list, gene_short_name, log2.transform=TRUE, min.pct.cells=10, downsample=FALSE, seed=1) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df.gene.norm <- MarvelObject$gene.norm.matrix
    cell.group.list <- cell.group.list
    gene_short_name <- gene_short_name
    min.pct.cells <- min.pct.cells
    downsample <- downsample
    seed <- seed
    log2.transform <- log2.transform
    
    # Example arguments
    #MarvelObject <- marvel
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #cell.group.list <- cell.group.list
    #gene_short_name <- "TPM1"
    #min.pct.cells <- 10
    #downsample <- TRUE
    #seed <- 1
    
    ##########################################################################
    
    # Downsample
    .list <- list()
    
    set.seed(seed)
    
    if(downsample==TRUE){
        
        # Find lowest common demoninator
        n.cells <- min(sapply(cell.group.list, length))
        
        # Track progress
        message(paste("Downsampling to ", n.cells, " cells per group", sep=""))
        
        # Downsample
        for(i in 1:length(cell.group.list)) {
            
            cell.ids <- cell.group.list[[i]]
            
            cell.ids <- sample(cell.ids, size=n.cells, replace=FALSE)
            
            .list[[i]] <- cell.ids
            
        }
        
        names(.list) <- names(cell.group.list)
        cell.group.list <- .list
            
    }

    # Subset relevant gene
    df.gene.norm <- df.gene.norm[gene_short_name, ]
    df.gene.norm <- data.frame("cell.id"=names(df.gene.norm),
                               "exp.gene.norm"=as.numeric(df.gene.norm),
                               stringsAsFactors=FALSE
                               )
    row.names(df.gene.norm) <- NULL
    
    # Tranform values
    if(log2.transform==TRUE) {
        
        df.gene.norm$exp.gene.norm <- log2(df.gene.norm$exp.gene.norm + 1)
    
    }
    
    # Annotate cell groups
        # Create ref table
        .list <- list()
        
        for(i in 1:length(cell.group.list)) {
            
            . <- data.frame("cell.id"=cell.group.list[[i]],
                            "group"=names(cell.group.list)[i]
                            )
                            
            .list[[i]] <- .
            
        }
        
        ref <- do.call(rbind.data.frame, .list)
        
        # Annotate
        df.gene.norm <- join(df.gene.norm, ref, by="cell.id", type="left")
        df.gene.norm <- df.gene.norm[!is.na(df.gene.norm$group), ]
        
        # Set factor levels
        df.gene.norm$group <- factor(df.gene.norm$group, levels=names(cell.group.list))
        df.gene.norm <- df.gene.norm[order(df.gene.norm$group), ]
        
    # Compute mean expression, % cells expressed
    groups <- levels(df.gene.norm$group)
    
    mean.expr <- NULL
    pct.cells.expr <- NULL
    
    for(i in 1:length(groups)) {
        
        # Define cell group
        group <- groups[i]
        
        # Subset
        df.gene.norm.small <- df.gene.norm[which(df.gene.norm$group==group), ]
        
        # Compute mean
        mean.expr[i] <- mean(df.gene.norm.small$exp.gene.norm)
        
        # Compute % cells expressing gene
        pct.cells.expr[i] <- sum(df.gene.norm.small$exp.gene.norm != 0) / nrow(df.gene.norm.small) * 100
        
    }
    
    results <- data.frame("group"=groups,
                          "mean.expr"=mean.expr,
                          "pct.cells.expr"=pct.cells.expr,
                          stringsAsFactors=FALSE
                          )
    
    # Censor lowly expressing cell types
    results$pct.cells.expr[which(results$pct.cells.expr < min.pct.cells)] <- NA
     
    # Set factor levels
    results$group <- factor(results$group, levels=names(cell.group.list))
    results <- results[order(results$group), ]
    
    # Bubble plot
        # Definition
        data <- results
        x <- as.character(rep(1, times=nrow(data)))
        y <- rep(c(nrow(data):1))
        z1 <- data$pct.cells.expr
        z2 <- data$mean.expr
        maintitle <- gene_short_name
        xtitle <- ""
        ytitle <- ""
        legendtitle.size <- "% cells"
        legendtitle.color <- "mean[log2(expr)]"
        labels.y <- data$group
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, size=z1, color=z2)) +
            scale_color_gradientn(colors=c("gray","cyan","green","yellow","red")) +
            scale_y_continuous(breaks=y, labels=labels.y) +
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
                )  +
                scale_size_area(breaks=c(25,50,75,100), limits=c(1, 100))
 
    ##################################################################
 
    # Save into new slots
    MarvelObject$adhocGene$Expression$Gene$Table <- results
    MarvelObject$adhocGene$Expression$Gene$Plot <- plot
    MarvelObject$adhocGene$cell.group.list <- cell.group.list
    MarvelObject$adhocGene$gene_short_name <- gene_short_name
    
    # Return final object
    return(MarvelObject)
            
}


