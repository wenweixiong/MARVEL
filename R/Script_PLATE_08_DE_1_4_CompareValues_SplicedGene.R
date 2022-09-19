#' @title Differential gene expression analysis for differentially spliced genes
#'
#' @description Performs differential gene expression analysis between 2 groups of cells only on differentially spliced genes.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2.
#' @param downsample Logical value. If set to \code{TRUE}, the number of cells in each cell group will be downsampled to the sample size of the smaller cell group so that both cell groups will have the sample size prior to differential expression analysis. Default is \code{FALSE}.
#' @param psi.method Vector of character string(s). To include significant events from these method(s) for differential gene expression analysis.
#' @param psi.pval Vector of numeric value(s). The adjusted p-value, below which, the splicing event is considered differentially spliced, and the corresponding genes will be included for differential gene expression analysis.
#' @param psi.delta Numeric value. The absolute difference in mean PSI values between \code{cell.group.g1} and \code{cell.group.g1}, above which, the splicing event is considered differentially spliced, and the corresponding genes will be included for differential gene expression analysis.
#' @param method.de.gene Character string. Same as \code{method} in \code{CompareValues} function.
#' @param method.adjust.de.gene Character string. Same as \code{method} in \code{CompareValues} function.
#' @param show.progress Logical value. If set to \code{TRUE}, progress bar will be displayed so that users can estimate the time needed for differential analysis. Default value is \code{TRUE}.
# @param use.downsampled.data Logical value. If set to \code{TRUE}, downsampled cell groups based on number of genes detected will be used and the sample IDs specified in \code{cell.group.g1} and \code{cell.group.g2} options will be overriden. The function \code{DownsampleByGenes} would need to be executed first if this option is set to \code{TRUE}. Default value is \code{FALSE}.
#'
#' @return An object of class S3 new slot \code{MarvelObject$DE$Exp$Table}.
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import utils
#' @export

CompareValues.Exp.Spliced <- function(MarvelObject, cell.group.g1=NULL, cell.group.g2=NULL, psi.method, psi.pval, psi.delta, method.de.gene="wilcox", method.adjust.de.gene="fdr", downsample=FALSE, show.progress=TRUE, use.downsampled.data=FALSE) {

    # Define arguments
    MarvelObject <- MarvelObject
    cell.group.g1 <- cell.group.g1
    cell.group.g2 <- cell.group.g2
    psi.method <- psi.method
    psi.pval <- psi.pval
    psi.delta <- psi.delta
    method.de.gene <- method.de.gene
    method.adjust.de.gene <- method.adjust.de.gene
    downsample <- downsample
    show.progress <- show.progress
    use.downsampled.data <- use.downsampled.data
    
    # Example arguments
    #MarvelObject <- marvel
    #cell.group.g1 <- cell.group.g1
    #cell.group.g2 <- cell.group.g2
    #psi.method <- c("ad", "dts")
    #psi.pval <- c(0.1, 0.1)
    #psi.delta <- 10
    #method.de.gene <- "wilcox"
    #method.adjust.de.gene <- "fdr"
    #downsample <- FALSE
    #show.progress <- TRUE
    #use.downsampled.data <- TRUE
    
    ##################################################
    
    # Tabulate sig events
    .list <- list()
    
    for(i in 1:length(psi.method)) {
    
        # Subset relevent splicing DE results
        de.psi <- MarvelObject$DE$PSI$Table[[psi.method[i]]]

        # Subset sig events
        index <- which(abs(de.psi$mean.diff) > psi.delta & de.psi$p.val.adj < psi.pval[i] & de.psi$outlier==FALSE)
        de.psi <- de.psi[index, ]
        
        # Subset gene metadata
        cols <- c("gene_id", "gene_short_name", "gene_type")
        de.psi <- de.psi[, cols]
        
        # Save into list
        .list[[i]] <- de.psi
        
    }
    
    df <- do.call(rbind.data.frame, .list)
    df <- unique(df)
    
    # Use downsampled cells based on no. of genes expressed
    if(use.downsampled.data==TRUE){
        
        df.downsampled <- MarvelObject$DE$Exp$Downsampled.Data$Table
        
        cell.group.g1 <- df.downsampled[which(df.downsampled$cell.group=="cell.group.g1"), "sample.id"]
        
        cell.group.g2 <- df.downsampled[which(df.downsampled$cell.group=="cell.group.g2"), "sample.id"]
        
    }
    
    # Perform DE analysis
    object <- CompareValues(MarvelObject=marvel,
                            cell.group.g1=cell.group.g1,
                            cell.group.g2=cell.group.g2,
                            min.cells=3,
                            method=method.de.gene,
                            method.adjust=method.adjust.de.gene,
                            level="gene",
                            custom.gene_ids=df$gene_id,
                            downsample=downsample,
                            show.progress=show.progress
                            )
                            
    ############################################################
    
    # Save to new slot
    MarvelObject$DE$Exp.Spliced$Table <- object$DE$Exp.Custom$Table
    
    # Return final object
    return(MarvelObject)
        
}
