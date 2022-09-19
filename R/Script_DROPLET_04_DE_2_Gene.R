#' @title Differential gene expression analysis
#'
#' @description Performs differential gene expression analysis between two groups of cells. Only among cells and genes previously included for splice junction analysis.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.SJ.10x} function.
#' @param log2.transform Logical value. If set to \code{TRUE} (default), normalised gene expression values will be off-set by 1 and then log2-transformed prior to analysis.
#' @param show.progress Logical value. If set to \code{TRUE} (default), the progress bar will appear.
#'
#' @return An object of class S3 with a updated slot \code{MarvelObject$DE$SJ$Table}.
#'
#' @importFrom plyr join
#' @import utils
#'
#' @export

CompareValues.Genes.10x <- function(MarvelObject, log2.transform=TRUE, show.progress=TRUE) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    de <- MarvelObject$DE$SJ$Table
    df.gene.norm <- MarvelObject$gene.norm.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.group.g1 <- MarvelObject$DE$SJ$cell.group.g1
    cell.group.g2 <- MarvelObject$DE$SJ$cell.group.g2
    show.progress <- show.progress
    log2.transform <- log2.transform
    
    # Example arguments
    #MarvelObject <- marvel
    #de <- MarvelObject$DE$SJ$Table
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #cell.group.g1 <- MarvelObject$DE$SJ$cell.group.g1
    #cell.group.g2 <- MarvelObject$DE$SJ$cell.group.g2
    
    #################################################################
    
    # Subset genes of SJ included for DE
    gene_short_names <- unique(de$gene_short_name)
    
    # Subset cell group 1
    df.gene.norm.g1 <- df.gene.norm[gene_short_names, cell.group.g1]
    
    # Subset cell group 2
    df.gene.norm.g2 <- df.gene.norm[gene_short_names, cell.group.g2]
    
    # log2 transform values
    if(log2.transform==TRUE) {
        
        df.gene.norm.g1 <- log2(df.gene.norm.g1 + 1)
        df.gene.norm.g2 <- log2(df.gene.norm.g2 + 1)
            
    }
    
    # Tabulate expression: Group 1
        # Compute group size
        n.cells.total <- length(cell.group.g1)
        
        # Compute % expr cells
        n.cells.expr.gene.norm <- apply(df.gene.norm.g1, 1, function(x) {sum(x != 0)})
        pct.cells.expr.gene.norm <- round(n.cells.expr.gene.norm/n.cells.total * 100, digits=2)
        
        # Compute mean expr
        mean.expr.gene.norm <- apply(df.gene.norm.g1, 1, function(x) {mean(x)})
        
        # Tabulate results
        results <- data.frame("gene_short_name"=gene_short_names,
                              "n.cells.total.norm"=n.cells.total,
                              "n.cells.expr.gene.norm"=n.cells.expr.gene.norm,
                              "pct.cells.expr.gene.norm"=pct.cells.expr.gene.norm,
                              mean.expr.gene.norm=mean.expr.gene.norm,
                              stringsAsFactors=FALSE
                              )
        row.names(results) <- NULL
        
        # Indicate group
        names(results)[which(names(results) != "gene_short_name")] <- paste(names(results)[which(names(results) != "gene_short_name")], ".g1", sep="")
        
        # Save as new object
        results.g1 <- results
        
    # Tabulate expression: Group 1
        # Compute group size
        n.cells.total <- length(cell.group.g2)
        
        # Compute % expr cells
        n.cells.expr.gene.norm <- apply(df.gene.norm.g2, 1, function(x) {sum(x != 0)})
        pct.cells.expr.gene.norm <- round(n.cells.expr.gene.norm/n.cells.total * 100, digits=2)
        
        # Compute mean expr
        mean.expr.gene.norm <- apply(df.gene.norm.g2, 1, function(x) {mean(x)})
        
        # Tabulate results
        results <- data.frame("gene_short_name"=gene_short_names,
                              "n.cells.total.norm"=n.cells.total,
                              "n.cells.expr.gene.norm"=n.cells.expr.gene.norm,
                              "pct.cells.expr.gene.norm"=pct.cells.expr.gene.norm,
                              mean.expr.gene.norm=mean.expr.gene.norm,
                              stringsAsFactors=FALSE
                              )
        row.names(results) <- NULL
        
        # Indicate group
        names(results)[which(names(results) != "gene_short_name")] <- paste(names(results)[which(names(results) != "gene_short_name")], ".g2", sep="")
        
        # Save as new object
        results.g2 <- results
    
    # Merge group 1, 2
    index.l <- table(results.g1$gene_short_name==results.g2$gene_short_name)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
    
        results.g2$gene_short_name <- NULL
        results <- cbind.data.frame(results.g1, results.g2)
        
    } else {
        
        return(print("Error in merging tables from Group 1 and 2"))

    }
    
    # Compute log2fc
    results$log2fc.gene.norm <- results$mean.expr.gene.norm.g2 - results$mean.expr.gene.norm.g1
    
    # Compute p-values
    print("Performing Wilcox rank sum test...")
    
    if(show.progress==TRUE) {
        
        pb <- txtProgressBar(1, length(gene_short_names), style=3)
    
    }
    
    pval <- NULL
    
    for(i in 1:length(gene_short_names)) {
        
        # Retrieve values
        values.g1 <- as.numeric(df.gene.norm.g1[gene_short_names[i], ])
        values.g2 <- as.numeric(df.gene.norm.g2[gene_short_names[i], ])
        
        # Wilcox
        pval[i] <- wilcox.test(values.g1, values.g2)$p.value
        
        if(show.progress==TRUE) {
            
            # Track progress
            setTxtProgressBar(pb, i)
            
        }
        
    }
    
    results$pval.gene.norm <- pval
    
    # Adjust for multiple testing
    results$pval.adj.gene.norm <- p.adjust(results$pval.gene.norm, method="fdr", n=length(results$pval.gene.norm))
    
    ################################################################
 
    # Update SJ DE table
    de <- join(de, results, by="gene_short_name", type="left")
 
    ################################################################
    # Save into new slot
    MarvelObject$DE$SJ$Table <- de
    
    # Return final object
    return(MarvelObject)
            
}


