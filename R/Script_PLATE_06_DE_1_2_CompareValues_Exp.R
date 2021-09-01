#' @title Differential Gene Expression Analysis
#'
#' @description
#' \code{CompareValues.Exp} performs differentially gene expression analysis between 2 groups of cells.
#'
#' @details
#' This function compares the gene expression values between 2 groups of cells.
#'
#' @param MarvelObject S3 object generated from \code{ComputePSI} function.
#' @param cell.type.columns.1 Character string. To indicate which columns in the \code{MarvelObject$GenePheno} slot to refer to when filtering samples (cells) for analysis. This is for the first group of cells (reference group).
#' @param cell.type.variables.1 List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{cell.type.columns} argument. This is for the first group of cells (reference group).
#' @param cell.type.columns.2 Character string. To indicate which columns in the \code{MarvelObject$GenePheno} slot to refer to when filtering samples (cells) for analysis. This is for the 2nd group of cells (non-reference group).
#' @param cell.type.variables.2 List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{cell.type.columns} argument. This is for the 2nd group of cells  (non-reference group).
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event for the event to be included for differential splicing analysis.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"kuiper"}, \code{"ad"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Kuiper, Anderson-Darling, Wilcox, and t-test, respectively.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$DE$Gene}
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- CompareValues.Exp(MarvelObject=marvel,
#'                             cell.type.columns.1=c("cell.type"),
#'                             cell.type.variables.1=list("iPSC"),
#'                             cell.type.columns.2=c("cell.type"),
#'                             cell.type.variables.2=list("Endoderm"),
#'                             n.cells=2,
#'                             method="t.test",
#'                             method.adjust="fdr"
#'                             )
#'
#' # Check output
#' marvel$DE$Exp$Table[1:5, ]

CompareValues.Exp <- function(MarvelObject, cell.type.columns.1, cell.type.variables.1, cell.type.columns.2, cell.type.variables.2, n.cells, method, method.adjust) {

    # Define arguments
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$GenePheno
    df.feature <- MarvelObject$GeneFeature
    cell.type.columns.1 <- cell.type.columns.1
    cell.type.variables.1 <- cell.type.variables.1
    cell.type.columns.2 <- cell.type.columns.2
    cell.type.variables.2 <- cell.type.variables.2
    n.cells <- n.cells
    method <- method
    method.adjust <- method.adjust
    
    # Define arguments
    #df <- marvel$Exp
    #df.pheno <- marvel$GenePheno
    #df.feature <- marvel$GeneFeature
    #cell.type.columns.1 <- c("cell.type")
    #cell.type.variables.1 <- list("iPSC")
    #cell.type.columns.2 <- c("cell.type")
    #cell.type.variables.2 <-list("Endoderm")
    #n.cells <- 3
    #method <- "wilcox"
    #method.adjust <- "fdr"
        
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL

    # Subset events with sufficient cells
        # Group 1
        .list <- list()
        
        for(i in 1:length(cell.type.columns.1)) {
            
            .list[[i]] <- df.pheno[which(df.pheno[[cell.type.columns.1[i]]] %in% cell.type.variables.1[[i]]), "sample.id"]
            
        }
        
        sample.ids.1 <- Reduce(intersect, .list)
        
        df.small <- df[, which(names(df) %in% sample.ids.1)]
        . <- apply(df.small, 1, function(x) {sum(x > 0)})
        gene_ids.1 <- names(.)[which(. > n.cells)]
        
        # Group 2
        .list <- list()
        
        for(i in 1:length(cell.type.columns.2)) {
            
            .list[[i]] <- df.pheno[which(df.pheno[[cell.type.columns.2[i]]] %in% cell.type.variables.2[[i]]), "sample.id"]
            
        }
        
        sample.ids.2 <- Reduce(intersect, .list)
        
        df.small <- df[, which(names(df) %in% sample.ids.2)]
        . <- apply(df.small, 1, function(x) {sum(x > 0)})
        gene_ids.2 <- names(.)[which(. > n.cells)]
       
        # Subset overlaps
        #overlap <- intersect(gene_ids.1, gene_ids.2)
        #df.feature <- df.feature[which(df.feature$gene_id %in% overlap), ]
        #df <- df[overlap, ]
        
        # Subset union
        all <- unique(c(gene_ids.1, gene_ids.2))
        df.feature <- df.feature[which(df.feature$gene_id %in% all), ]
        df <- df[df.feature$gene_id, ]
        
    # Statistical test
    gene_ids <- df.feature$gene_id
    
    n.cells.x <- NULL
    n.cells.y <- NULL
    mean.x <- NULL
    mean.y <- NULL
    log2fc <- NULL
    p.val <- NULL

    #pb <- txtProgressBar(1, length(gene_ids), style=3)

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
        if(method=="ks") {
        
            p.val[i] <- ks.test(x, y)$p.value
        
        } else if(method=="wilcox") {
        
            p.val[i] <- wilcox.test(x, y)$p.value
        
         
        } else {
        
            p.val[i] <- t.test(x, y)$p.value
        
        }
        
        # Track progress
        #setTxtProgressBar(pb, i)
        
    }
        
    # Save into data frame
    results <- data.frame("gene_id"=gene_ids,
                          "n.cells.g1"=n.cells.x, "n.cells.g2"=n.cells.y,
                          "mean.g1"=mean.x, "mean.g2"=mean.y,
                          "log2fc"=log2fc,
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
    print(paste(sum(results$p.val.adj < 0.10), " DE genes < 0.10 adjusted p-value", sep=""))
    print(paste(sum(results$p.val.adj < 0.05), " DE genes < 0.05 adjusted p-value", sep=""))
    print(paste(sum(results$p.val.adj < 0.01), " DE genes < 0.01 adjusted p-value", sep=""))
    
    # Save to new slot
    MarvelObject$DE$Exp$Table <- results
  
    return(MarvelObject)
        
}
