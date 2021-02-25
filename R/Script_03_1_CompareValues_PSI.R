#' @title Differential splicing analysis
#'
#' @description
#' \code{CompareValues.PSI} performs differentially splicing analysis between 2 groups of cells.
#'
#' @details
#' This function compares the percent spliced-in (PSI) values between 2 groups of cells.
#'
#' @param MarvelObject S3 object generated from \code{ComputePSI} function.
#' @param cell.types Character string. To indicate which 2 groups of cells that will be used for differential splicing analysis. Group names should match those in \code{cell.type} column of \code{$SplicePheno} slot.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event for the event to be included for differential splicing analysis.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Wilcox, and t-test, respectively. We advice \code{"ks"} for PSI comparison while \code{"wilcox"} or \code{"t.test"} for gene expression comparison.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{$DE$PSI}
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' marvel <- CompareValues.PSI(MarvelObject=marvel,
#'                             cell.types=c("iPSC", "Endoderm"),
#'                             n.cells=25,
#'                             method="t.test",
#'                             method.adjust="fdr"
#'                             )
#'
#' marvel$DE$PSI[1:5, ]

CompareValues.PSI <- function(MarvelObject, cell.types, n.cells, method, method.adjust) {

    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.types <- cell.types
    n.cells <- n.cells
    method <- method
    method.adjust <- method.adjust
    
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Subset overlapping samples in matrix and pheno file
    df <- df[, which(names(df) %in% df.pheno$sample.id)]

    # Subset events with sufficient cells
        # Group 1
        sample.ids <- df.pheno[which(df.pheno$cell.type %in% cell.types[1]), "sample.id"]
        df.small <- df[, which(names(df) %in% sample.ids)]
        . <- apply(df.small, 1, function(x) {sum(!is.na(x))})
        tran_ids.1 <- names(.)[which(. >= n.cells)]
        
        # Group 2
        sample.ids <- df.pheno[which(df.pheno$cell.type %in% cell.types[2]), "sample.id"]
        df.small <- df[, which(names(df) %in% sample.ids)]
        . <- apply(df.small, 1, function(x) {sum(!is.na(x))})
        tran_ids.2 <- names(.)[which(. >= n.cells)]
       
        # Subset overlaps
        overlap <- intersect(tran_ids.1, tran_ids.2)
        df.feature <- df.feature[which(df.feature$tran_id %in% overlap), ]
        df <- df[overlap, ]
        
    # Check if matrix column and rows align with metadata
        # Column
        index.check <- which(unique((names(df)==df.pheno$sample.id))==FALSE)
        
        if(length(index.check)==0) {
            
            print("Checking... Matrix column (sample) names match sample metadata")
            
        } else {
            
            print("Checking... Matrix column (sample) names DO NOT match sample metadata")
            
        }
        
        # Row
        index.check <- which(unique((row.names(df)==df.feature$tran.id))==FALSE)
        
        if(length(index.check)==0) {
            
            print("Checking... Matrix row (feature) names match feature metadata")
            
        } else {
            
            print("Checking... Matrix row (feature) names DO NOT match feature metadata")
            
        }
        
    # Statistical test
    tran_ids <- row.names(df)
    
    n.cells.x <- NULL
    n.cells.y <- NULL
    mean.x <- NULL
    mean.y <- NULL
    mean.diff <- NULL
    p.val <- NULL

    #pb <- txtProgressBar(1, length(tran_ids), style=3)

    for(i in 1:length(tran_ids)) {

        # Subset event
        . <- df[tran_ids[i], ]
        . <- as.data.frame(t(.))
        . <- na.omit(.)
        names(.) <- "psi"
        .$sample.id <- row.names(.)
        row.names(.) <- NULL
        
        # Subset relevant cell types
        . <- join(., df.pheno[,c("sample.id", "cell.type")], by="sample.id", type="left")
        . <- .[which(.$cell.type %in% c(cell.types[1], cell.types[2])), ]
        
        # Retrieve values
        x <- .[which(.$cell.type==cell.types[1]), "psi"]
        y <- .[which(.$cell.type==cell.types[2]), "psi"]
        
        # Compute statistics
        n.cells.x[i] <- length(x)
        n.cells.y[i] <- length(y)
        mean.x[i] <- mean(x)
        mean.y[i] <- mean(y)
        mean.diff[i] <- mean(y) - mean(x)
        
        # Statistical test
        if(method=="ks") {
        
            p.val[i] <- ks.test(x, y)$p.value
        
        } else if(method=="wilcox") {
        
            p.val[i] <- wilcox.test(psi ~ cell.type, .)$p.value

         
        } else {
        
            p.val[i] <- t.test(psi ~ cell.type, .)$p.value

        }
        
        # Track progress
        #setTxtProgressBar(pb, i)
        
    }
        
    # Save into data frame
    results <- data.frame("tran_id"=tran_ids,
                          "n.cells.g1"=n.cells.x, "n.cells.g2"=n.cells.y,
                          "mean.g1"=mean.x, "mean.g2"=mean.y,
                          "mean.diff"=mean.diff,
                          "p.val"=p.val,
                          stringsAsFactors=FALSE)
    
    # Reorder by p-value
    results <- results[order(results$p.val), ]
    
    # Adjust for multiple testing
    results$p.val.adj <- p.adjust(results$p.val, method=method.adjust, n=length(results$p.val))
    
    # Annotate with feature metadata
    results <- join(results, df.feature, by="tran_id", type="left")
    cols.1 <- names(df.feature)
    cols.2 <- setdiff(names(results), names(df.feature))
    results <- results[,c(cols.1, cols.2)]
                          
    # Save to new slot
    MarvelObject$DE$PSI <- results
  
    return(MarvelObject)
        
}
