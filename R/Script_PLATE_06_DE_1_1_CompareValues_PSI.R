#' @title Differential Splicing Analysis
#'
#' @description
#' \code{CompareValues.PSI} performs differentially splicing analysis between 2 groups of cells.
#'
#' @details
#' This function compares the percent spliced-in (PSI) values between 2 groups of cells.
#'
#' @param MarvelObject S3 object generated from \code{ComputePSI} function.
#' @param cell.type.columns.1 Character string. To indicate which columns in the \code{MarvelObject$SplicePheno} slot to refer to when filtering samples (cells) for analysis. This is for the first group of cells (reference group).
#' @param cell.type.variables.1 List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{cell.type.columns.1} argument. This is for the first group of cells (reference group).
#' @param cell.type.columns.2 Character string. To indicate which columns in the \code{MarvelObject$SplicePheno} slot to refer to when filtering samples (cells) for analysis. This is for the 2nd group of cells (non-reference group).
#' @param cell.type.variables.2 List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{cell.type.columns.2} argument. This is for the 2nd group of cells  (non-reference group).
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event for the event to be included for differential splicing analysis.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"kuiper"}, \code{"ad"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Kuiper, Anderson-Darling, Wilcox, and t-test, respectively.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param event.type Character string. Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$DE$PSI}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import kuiper.2samp
#' @import kSamples
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- CompareValues.PSI(MarvelObject=marvel,
#'                             cell.type.columns.1=c("cell.type"),
#'                             cell.type.variables.1=list("iPSC"),
#'                             cell.type.columns.2=c("cell.type"),
#'                             cell.type.variables.2=list("Endoderm"),
#'                             n.cells=2,
#'                             method="t.test",
#'                             method.adjust="fdr",
#'                             event.type=c("SE", "MXE", "RI", "A5SS", "A3SS")
#'                             )
#'
#' # Check output
#' marvel$DE$PSI$Table[1:5, ]


CompareValues.PSI <- function(MarvelObject, cell.type.columns.1, cell.type.variables.1, cell.type.columns.2, cell.type.variables.2, n.cells, method, method.adjust, event.type) {

    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.type.columns.1 <- cell.type.columns.1
    cell.type.variables.1 <- cell.type.variables.1
    cell.type.columns.2 <- cell.type.columns.2
    cell.type.variables.2 <- cell.type.variables.2
    n.cells <- n.cells
    method <- method
    method.adjust <- method.adjust
    event.type <- event.type
    
    # Example arguments
    #df <- do.call(rbind.data.frame, marvel$PSI)
    #df.pheno <- marvel$SplicePheno
    #df.feature <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
    #cell.type.columns.1 <- c("cell.type")
    #cell.type.variables.1 <- list(c("iPSC"))
    #cell.type.columns.2 <- c("cell.type")
    #cell.type.variables.2 <- list("Endoderm")
    #n.cells <- 2
    #method <- "t.test"
    #method.adjust <- "fdr"
    #event.type <- c("SE", "MXE", "RI", "A5SS", "A3SS")
    
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Subset relevant events
    df.feature <- df.feature[which(df.feature$event_type %in% event.type), ]
    df <- df[df.feature$tran_id,]
        
    # Subset sample IDs (for DE later) and events with sufficient cells
        # Group 1
        .list <- list()
        
        for(i in 1:length(cell.type.columns.1)) {
            
            .list[[i]] <- df.pheno[which(df.pheno[[cell.type.columns.1[i]]] %in% cell.type.variables.1[[i]]), "sample.id"]
            
        }
        
        sample.ids.1 <- Reduce(intersect, .list)
        
        df.small <- df[, which(names(df) %in% sample.ids.1)]
        . <- apply(df.small, 1, function(x) {sum(!is.na(x))})
        tran_ids.1 <- names(.)[which(. >= n.cells)]
        
        # Group 2
        .list <- list()
        
        for(i in 1:length(cell.type.columns.2)) {
            
            .list[[i]] <- df.pheno[which(df.pheno[[cell.type.columns.2[i]]] %in% cell.type.variables.2[[i]]), "sample.id"]
            
        }
        
        sample.ids.2 <- Reduce(intersect, .list)
        
        df.small <- df[, which(names(df) %in% sample.ids.2)]
        . <- apply(df.small, 1, function(x) {sum(!is.na(x))})
        tran_ids.2 <- names(.)[which(. >= n.cells)]
       
        # Subset overlaps
        overlap <- intersect(tran_ids.1, tran_ids.2)
        df.feature <- df.feature[which(df.feature$tran_id %in% overlap), ]
        df <- df[overlap, ]
                        
    # Statistical test
    tran_ids <- row.names(df)
    
    n.cells.x <- NULL
    n.cells.y <- NULL
    mean.x <- NULL
    mean.y <- NULL
    mean.diff <- NULL
    p.val <- NULL
    statistic <- NULL
    
    #pb <- txtProgressBar(1, length(tran_ids), style=3)

    for(i in 1:length(tran_ids)) {

        # Subset event
        . <- df[tran_ids[i], ]
        . <- as.data.frame(t(.))
        . <- na.omit(.)
        names(.) <- "psi"
        .$sample.id <- row.names(.)
        row.names(.) <- NULL
        
        # Retrieve values
        x <- .[which(.$sample.id %in% sample.ids.1), "psi"]
        y <- .[which(.$sample.id %in% sample.ids.2), "psi"]
        
        # Compute statistics
        n.cells.x[i] <- length(x)
        n.cells.y[i] <- length(y)
        mean.x[i] <- mean(x)
        mean.y[i] <- mean(y)
        mean.diff[i] <- mean(y) - mean(x)
        
        # Statistical test
        if(method=="wilcox") {
            
            statistic[i] <- NA
            p.val[i] <- wilcox.test(x, y)$p.value

         
        } else if(method=="t.test"){
        
            statistic <- t.test(x, y)$statistic
            p.val[i] <- t.test(x, y)$p.value

        } else if(method=="ks") {
            
            statistic[i] <- ks.test(x, y)$statistic
            p.val[i] <- ks.test(x, y)$p.value
            
        } else if(method=="kuiper") {
            
            statistic[i] <- kuiper.2samp(x, y)$Kuiper.statistic
            p.val[i] <- kuiper.2samp(x, y)$p.value
            
        } else if(method=="ad") {
            
            error.check <- tryCatch(ad.test(x, y), error=function(err) "Error")
            
            if(error.check == "Error") {
                
                statistic[i] <- 0
                p.val[i] <- 1
                
            } else {
                
                statistic[i] <- ad.test(x, y, method="asymptotic")$ad[1,1]
                p.val[i] <- ad.test(x, y, method="asymptotic")$ad[1,3]
                
            }
            
        }
        
        # Track progress
        #setTxtProgressBar(pb, i)
        #print(i)
        
    }
        
    # Save into data frame
    results <- data.frame("tran_id"=tran_ids,
                          "n.cells.g1"=n.cells.x, "n.cells.g2"=n.cells.y,
                          "mean.g1"=mean.x, "mean.g2"=mean.y,
                          "mean.diff"=mean.diff,
                          "statistic"=as.numeric(statistic),
                          "p.val"=p.val,
                          stringsAsFactors=FALSE
                          )
    
    # Reorder by p-value
    results <- results[which(!is.na(results$p.val)), ]
    results$statistic[which(results$p.val < 0)]  <- results$statistic[which(results$p.val < 0)] * -1
    results$p.val <- abs(results$p.val)
    results <- results[order(results$p.val), ]
    #results <- results[order(results$statistic, decreasing=TRUE), ]
    
    # Adjust for multiple testing
    results$p.val.adj <- p.adjust(results$p.val, method=method.adjust, n=length(results$p.val))
    
    # Annotate with feature metadata
    results <- join(results, df.feature, by="tran_id", type="left")
    cols.1 <- names(df.feature)
    cols.2 <- setdiff(names(results), names(df.feature))
    results <- results[,c(cols.1, cols.2)]
    
    # Report result summary
    print(paste(sum(results$p.val.adj < 0.10), " DE splicing events < 0.10 adjusted p-value", sep=""))
    print(paste(sum(results$p.val.adj < 0.05), " DE splicing events < 0.05 adjusted p-value", sep=""))
    print(paste(sum(results$p.val.adj < 0.01), " DE splicing events < 0.01 adjusted p-value", sep=""))
                          
    # Save to new slot
    MarvelObject$DE$PSI$Table <- results
    MarvelObject$DE$PSI$cell.type.columns.1 <- cell.type.columns.1
    MarvelObject$DE$PSI$cell.type.variables.1 <- cell.type.variables.1
    MarvelObject$DE$PSI$cell.type.columns.2 <- cell.type.columns.2
    MarvelObject$DE$PSI$cell.type.variables.2 <- cell.type.variables.2
    MarvelObject$DE$PSI$n.cells <- n.cells
    return(MarvelObject)
        
}
