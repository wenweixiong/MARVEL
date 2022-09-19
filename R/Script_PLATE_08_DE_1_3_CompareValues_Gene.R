#' @title Differential gene expression analysis
#'
#' @description Performs differential gene expression analysis between 2 groups of cells.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2.
#' @param downsample Logical value. If set to \code{TRUE}, the number of cells in each cell group will be downsampled to the sample size of the smaller cell group so that both cell groups will have the sample size prior to differential expression analysis. Default is \code{FALSE}.
#' @param min.cells Numeric value. The minimum no. of cells expressing the gene for the gene to be included for differential splicing analysis.
#' @param pct.cells Numeric value. The minimum no. of cells expressing the gene for the gene to be included for differential splicing analysis. If \code{pct.cells} is specified, then \code{pct.cells} will be used as threshold instead of \code{min.cells}.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"kuiper"}, \code{"ad"}, \code{"dts"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Kuiper, Anderson-Darling, DTS, Wilcox, and t-test, respectively.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param show.progress Logical value. If set to \code{TRUE}, progress bar will be displayed so that users can estimate the time needed for differential analysis. Default value is \code{TRUE}.
#' @param nboots Numeric value. When \code{method} set to \code{"dts"}, the number of bootstrap iterations for computing the p-value.
#' @param custom.gene_ids Character string. Instead of specified the genes to include for DE analysis with \code{min.cells}, users may input a custom vector of gene IDs to include for DE analysis.
# @param use.downsampled.data Logical value. If set to \code{TRUE}, downsampled cell groups based on number of genes detected will be used and the sample IDs specified in \code{cell.group.g1} and \code{cell.group.g2} options will be overriden. The function \code{DownsampleByGenes} would need to be executed first if this option is set to \code{TRUE}. Default value is \code{FALSE}.
#'
#' @return An object of class S3 new slot \code{MarvelObject$DE$Exp$Table}.
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import utils
#' @export

CompareValues.Exp <- function(MarvelObject, cell.group.g1=NULL, cell.group.g2=NULL, downsample=FALSE, min.cells=25, pct.cells=NULL, method, method.adjust, show.progress=TRUE, nboots=1000, custom.gene_ids=NULL, use.downsampled.data=FALSE) {

    # Define arguments
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- MarvelObject$GeneFeature
    cell.group.g1 <- cell.group.g1
    cell.group.g2 <- cell.group.g2
    downsample <- downsample
    min.cells <- min.cells
    pct.cells <- pct.cells
    method <- method
    method.adjust <- method.adjust
    show.progress <- show.progress
    nboots <- nboots
    custom.gene_ids <- custom.gene_ids
    use.downsampled.data <- use.downsampled.data
    
    # Define arguments
    #df <- marvel$Exp
    #df.pheno <- marvel$SplicePheno
    #df.feature <- marvel$GeneFeature
    #cell.group.g1 <- cell.group.g1
    #cell.group.g2 <- cell.group.g2
    #downsample <- TRUE
    #min.cells <- 25
    #pct.cells <- NULL
    #method <- "wilcox"
    #method.adjust <- "fdr"
    #show.progress <- TRUE
    #custom.gene_ids <- custom.gene_ids
    
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
  
    # Use downsampled cells based on no. of genes expressed
    if(use.downsampled.data==TRUE){
      
        df.downsampled <- MarvelObject$DE$Exp$Downsampled.Data$Table
      
        cell.group.g1 <- df.downsampled[which(df.downsampled$cell.group=="cell.group.g1"), "sample.id"]
      
        cell.group.g2 <- df.downsampled[which(df.downsampled$cell.group=="cell.group.g2"), "sample.id"]
      
    }
  
    # Retrieve sample ids
        # Group 1
        sample.ids.1 <- cell.group.g1
        
        # Group 2
        sample.ids.2 <- cell.group.g2
    
    # Downsample
    if(downsample==TRUE) {
        
        # Retrieve lowest denominator
        n.cells.downsample <- min(length(sample.ids.1), length(sample.ids.2))
        
        # Downsample
        set.seed(1)
        sample.ids.1.small <- sample(sample.ids.1, size=n.cells.downsample, replace=FALSE)
        sample.ids.2.small <- sample(sample.ids.2, size=n.cells.downsample, replace=FALSE)
        
        # Subset cells
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% c(sample.ids.1.small, sample.ids.2.small)), ]
        df <- df[, df.pheno$sample.id]
        
        # Track progress
        print(paste(length(sample.ids.1), " cells found in Group 1", sep=""))
        print(paste(length(sample.ids.2), " cells found in Group 2", sep=""))
        print(paste("Both Group 1 and 2 downsampled to ", n.cells.downsample, " cells", sep=""))
        
    }
    
    # Subset events with sufficient cells
    if(is.null(custom.gene_ids[1])) {
        
        # Group 1
        df.small <- df[, which(names(df) %in% sample.ids.1)]
        . <- apply(df.small, 1, function(x) {sum(x > 0)})
        
        if(is.null(pct.cells)) {
            
            gene_ids.1 <- names(.)[which(. >= min.cells)]
        
        } else {
            
            . <- ./nrow(df.pheno) * 100
            
            gene_ids.1 <- names(.)[which(. >= pct.cells)]
            
        }
        
        # Group 2
        df.small <- df[, which(names(df) %in% sample.ids.2)]
        . <- apply(df.small, 1, function(x) {sum(x > 0)})
        
        if(is.null(pct.cells)) {
            
            gene_ids.2 <- names(.)[which(. >= min.cells)]
        
        } else {
            
            . <- ./nrow(df.pheno) * 100
            
            gene_ids.2 <- names(.)[which(. >= pct.cells)]
            
        }
       
        # Subset overlaps
        #overlap <- intersect(gene_ids.1, gene_ids.2)
        #df.feature <- df.feature[which(df.feature$gene_id %in% overlap), ]
        #df <- df[overlap, ]
                
        # Subset union
        all <- unique(c(gene_ids.1, gene_ids.2))
        df.feature <- df.feature[which(df.feature$gene_id %in% all), ]
        df <- df[df.feature$gene_id, ]
        
        # Report progress
        print(paste(length(gene_ids.1), " expressed genes identified in Group 1", sep=""))
        print(paste(length(gene_ids.2), " expressed genes identified in Group 2", sep=""))
        print(paste(length(all), " expressed genes identified in EITHER Group 1 or Group 2", sep=""))
    
    } else {
        
        df.feature <- df.feature[which(df.feature$gene_id %in% custom.gene_ids), ]
        df <- df[df.feature$gene_id, ]
        print(paste(length(custom.gene_ids), " custom genes specified", sep=""))
        
    }
    
    # Statistical test
    gene_ids <- df.feature$gene_id
    
    n.cells.x <- NULL
    n.cells.y <- NULL
    mean.x <- NULL
    mean.y <- NULL
    log2fc <- NULL
    statistic <- NULL
    p.val <- NULL

    if(show.progress==TRUE) {
        
        pb <- txtProgressBar(1, length(gene_ids), style=3)
        
    }

    for(i in 1:length(gene_ids)) {

        # Subset event
        . <- df[gene_ids[i], ]
        . <- as.data.frame(t(.))
        . <- na.omit(.)
        names(.) <- "exp"
        .$sample.id <- row.names(.)
        row.names(.) <- NULL
                
        # Retrieve values
        x <- .[which(.$sample.id %in% sample.ids.1), "exp"]
        y <- .[which(.$sample.id %in% sample.ids.2), "exp"]
        
        # Compute statistics
        n.cells.x[i] <- length(which(x > 0))
        n.cells.y[i] <- length(which(y > 0))
        mean.x[i] <- mean(x)
        mean.y[i] <- mean(y)
        log2fc[i] <- mean(y) - mean(x)
        
        # Statistical test
        if(method=="wilcox") {
            
            statistic[i] <- NA
            p.val[i] <- wilcox.test(x, y)$p.value

         
        } else if(method=="t.test"){
        
            statistic[i] <- t.test(x, y)$statistic
            p.val[i] <- t.test(x, y)$p.value

        } else if(method=="ks") {
            
            statistic[i] <- ks.test(x, y)$statistic
            p.val[i] <- ks.test(x, y)$p.value
            
        } else if(method=="kuiper") {
            
            statistic[i] <- kuiper.2samp(x, y)$Kuiper.statistic
            p.val[i] <- kuiper.2samp(x, y)$p.value
            
        } else if(method=="ad") {
            
            error.check <- tryCatch(ad.test(x, y), error=function(err) "Error")
            
            if(error.check[1] == "Error") {
                
                statistic[i] <- 0
                p.val[i] <- 1
                
            } else {
                
                . <- ad.test(x, y, method="asymptotic")$ad
                
                statistic[i] <- .[1,1]
                p.val[i] <- .[1,3]
                
            }
            
        } else if(method=="dts"){
            
            . <- dts_test(x, y, nboots=nboots)
            
            statistic[i] <- .[1]
            p.val[i] <- .[2]
            
            
        }
        
        # Track progress
        if(show.progress==TRUE) {
            
            setTxtProgressBar(pb, i)
            
        }
        
    }
        
    # Save into data frame
    results <- data.frame("gene_id"=gene_ids,
                          "n.cells.g1"=n.cells.x, "n.cells.g2"=n.cells.y,
                          "mean.g1"=mean.x, "mean.g2"=mean.y,
                          "log2fc"=log2fc,
                          "statistic"=statistic,
                          "p.val"=p.val,
                          stringsAsFactors=FALSE)
    
    # Reorder by p-value
    results <- results[order(results$p.val), ]
    
    # Adjust for multiple testing
    results$p.val.adj <- p.adjust(results$p.val, method=method.adjust, n=length(results$p.val))
    
    # Annotate with feature metadata
    results <- join(results, df.feature, by="gene_id", type="left")
    cols.1 <- names(df.feature)
    cols.2 <- setdiff(names(results), names(df.feature))
    results <- results[,c(cols.1, cols.2)]
    
    # Report result summary
    #print(paste(sum(results$p.val.adj < 0.10), " DE genes < 0.10 adjusted p-value", sep=""))
    #print(paste(sum(results$p.val.adj < 0.05), " DE genes < 0.05 adjusted p-value", sep=""))
    #print(paste(sum(results$p.val.adj < 0.01), " DE genes < 0.01 adjusted p-value", sep=""))
    
    # Save to new slot
    if(is.null(custom.gene_ids[1])) {
        
        MarvelObject$DE$Exp$Table <- results
        
    } else {
        
        MarvelObject$DE$Exp.Custom$Table <- results
        
    }
  
    return(MarvelObject)
        
}
