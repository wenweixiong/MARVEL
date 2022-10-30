#' @title Differential gene expression analysis
#'
#' @description Performs differential gene expression analysis between two groups of cells. Only among cells and genes previously included for splice junction analysis.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.SJ.10x} function.
#' @param log2.transform Logical value. If set to \code{TRUE} (default), normalised gene expression values will be off-set by 1 and then log2-transformed prior to analysis. This option is automatically set to \code{TRUE} if \code{method} option is set to \code{"mast"}.
#' @param method Character string. Statistical test to compare the 2 groups of cells. Default is \code{"wilcox"} as recommended by Seurat. Another option is \code{"mast"}. If \code{"mast"} is specified, the log2fc and p-values will be corrected using the gene detection rate as per the \code{MAST} package tutorial.
#' @param show.progress Logical value. If set to \code{TRUE} (default), the progress bar will appear.
#' @param mast.method Character string. As per the \code{method} option of the \code{zlm} function from the \code{MAST} package. Default is \code{"bayesglm"}, other options are \code{"glm"} and \code{"glmer"}.
#' @param mast.ebayes Logical value. As per the \code{ebayes} option of the \code{zlm} function from the \code{MAST} package. Default is \code{TRUE}.
#'
#' @return An object of class S3 with a updated slot \code{MarvelObject$DE$SJ$Table}.
#'
#' @importFrom plyr join
#' @importFrom stats p.adjust p.adjust.methods
#' @importFrom utils txtProgressBar setTxtProgressBar
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
#' marvel.demo.10x <- CompareValues.Genes.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         show.progress=FALSE
#'                         )
#'
#' # Check output
#' head(marvel.demo.10x$DE$SJ$Table)

CompareValues.Genes.10x <- function(MarvelObject, log2.transform=TRUE, show.progress=TRUE, method="wilcox", mast.method="bayesglm", mast.ebayes=TRUE) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    de <- MarvelObject$DE$SJ$Table
    df.gene.norm <- MarvelObject$gene.norm.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.group.g1 <- MarvelObject$DE$SJ$cell.group.g1
    cell.group.g2 <- MarvelObject$DE$SJ$cell.group.g2
    
    # Example arguments
    #MarvelObject <- marvel
    #de <- MarvelObject$DE$SJ$Table
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #cell.group.g1 <- MarvelObject$DE$SJ$cell.group.g1
    #cell.group.g2 <- MarvelObject$DE$SJ$cell.group.g2
    #log2.transform <- TRUE
    #method <- "mast"
    #mast.method <- "bayesglm"
    #mast.ebayes <- TRUE
    
    #################################################################
    
    # Subset genes of SJ included for DE
    gene_short_names <- unique(de$gene_short_name)
    
    # Subset cell group 1
    df.gene.norm.g1 <- df.gene.norm[gene_short_names, cell.group.g1]
    
    # Subset cell group 2
    df.gene.norm.g2 <- df.gene.norm[gene_short_names, cell.group.g2]
    
    # log2 transform values
    if(log2.transform==TRUE | method=="mast") {
        
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
        
    # Tabulate expression: Group 2
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
        
        return(message("Error in merging tables from Group 1 and 2"))

    }
    
    #################################################################
    
    if(method=="wilcox"){
        
        # Compute log2fc
        results$log2fc.gene.norm <- results$mean.expr.gene.norm.g2 - results$mean.expr.gene.norm.g1
        
        # Compute p-values
        message("Performing Wilcox rank sum test...")
        
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
        
    } else if(method=="mast"){
        
        # Prepare cdata
            # Indicate group 1,2
            cdata <- data.frame("cell.id"=c(cell.group.g1, cell.group.g2))
            cdata$condition <- ifelse(cdata$cell.id %in% cell.group.g1, "g1", "g2")
            cdata$condition <- factor(cdata$condition, levels=c("g1", "g2"))
            
            # Retrieve expression values
            df.exp.master <- cbind(df.gene.norm.g1, df.gene.norm.g2)
            df.exp.master <- df.exp.master[,cdata$cell.id]
            
            # Compute n expressed genes
            . <- apply(df.exp.master, 2, function(x) {sum(x != 0)})
            . <- data.frame("cell.id"=names(.),
                            "n.genes"=as.numeric(.),
                            stringsAsFactors=FALSE
                            )
                            
            # Compute gene detection rate
            .$cngeneson <- scale(.$n.genes)
            cdata <- join(cdata, ., by="cell.id", type="left")
            
            # Format for MAST
            names(cdata)[which(names(cdata)=="cell.id")] <- "wellKey"
            row.names(cdata) <- cdata$wellKey
            
        # Prepare fdata
        fdata <- data.frame("primerid"=row.names(df.exp.master),
                            stringsAsFactors=FALSE
                            )
        row.names(fdata) <- fdata$primerid
        
        # Prepare SingleCellAssay object
        sca <- MAST::FromMatrix(as.matrix(df.exp.master), cdata, fdata)
        
        # Build model
        zlmCond <- MAST::zlm(~condition + cngeneson, sca, method=mast.method, ebayes=mast.ebayes, silent=TRUE)

        # Only test the condition coefficient
        summaryCond <- summary(zlmCond, doLRT='conditiong2')

        # Retrieve DE table
        summaryDt <- summaryCond$datatable
        #fcHurdle <- merge(summaryDt[contrast=='conditiong2' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditiong2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
        #fcHurdle <- as.data.frame(fcHurdle)
        
        index <- which(summaryDt$contrast=='conditiong2' & summaryDt$component=='H')
        summaryDt.small.1 <- summaryDt[index, c("primerid", "Pr(>Chisq)")]
        
        index <- which(summaryDt$contrast=='conditiong2' & summaryDt$component=='logFC')
        summaryDt.small.2 <- summaryDt[index, c("primerid", "coef", "ci.hi", "ci.lo")]
        
        fcHurdle <- join(summaryDt.small.1, summaryDt.small.2, by="primerid")
        
        fcHurdle <- as.data.frame(fcHurdle)
        
        # Subset relevant columns
        results. <- fcHurdle
        results. <- results.[,c("primerid", "Pr(>Chisq)", "coef")]
        names(results.) <- c("gene_short_name", "pval.gene.norm", "log2fc.gene.norm")
        results. <- results.[,c("gene_short_name", "log2fc.gene.norm", "pval.gene.norm")]
        
        # Annotate original result table
        results <- join(results, results., by="gene_short_name", type="left")
        
    }
    
    #################################################################

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


