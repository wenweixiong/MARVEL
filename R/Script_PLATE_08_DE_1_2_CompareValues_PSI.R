#' @title Differential splicing analysis
#'
#' @description Performs differentially splicing analysis between 2 groups of cells.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2.
#' @param downsample Logical value. If set to \code{TRUE}, the number of cells in each cell group will be downsampled to the sample size of the smaller cell group so that both cell groups will have the sample size prior to differential expression analysis. Default is \code{FALSE}.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event for the event to be included for differential splicing analysis.
#' @param pct.cells Numeric value. The minimum percentage of cells expressing the splicing event for the event to be included for differential splicing analysis. If \code{pct.cells} is specified, then \code{pct.cells} will be used as threshold instead of \code{min.cells}.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"kuiper"}, \code{"ad"}, \code{"dts"}, \code{"wilcox"}, \code{"t.test"}, and \code{"permutation"} for Kolmogorov-Smirnov, Kuiper, Anderson-Darling, DTS, Wilcox, t-test, and, permutation approach respectively.
#' @param n.permutations Numeric value. When \code{method} set to \code{"permutation"}, this argument indicates the number of permutations to perform for generating the null distribution for subsequent p-value inference. Default is \code{1000} times.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param event.type Character string. Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param show.progress Logical value. If set to \code{TRUE}, progress bar will be displayed so that users can estimate the time needed for differential analysis. Default value is \code{TRUE}.
#' @param nboots Numeric value. When \code{method} set to \code{"dts"}, the number of bootstrap iterations for computing the p-value.
#' @param annotate.outliers Numeric value. When set to \code{TRUE}, statistical difference in PSI values between the two cell groups that is driven by outlier cells will be annotated.
#' @param n.cells.outliers Numeric value. When \code{annotate.outliers} set to \code{TRUE}, the minimum number of cells with non-1 or non-0 PSI values for included-to-included or excluded-to-excluded modality change, respectively. The p-values will be re-coded to 1 when both cell groups have less than this minimum number of cells. This is to avoid false positive results.
#' @param assign.modality Logical value. If set to \code{TRUE} (default), modalities will be assigned to each cell group.
#' @param seed Numeric value. The seed number for the random number generator to ensure reproducibility during during down-sampling of cells when \code{downsample} set to \code{TRUE}, during permutation testing when \code{method} set to \code{"permutation"}, and during modality assignment which will be performed automatically.
#'
#' @return An object of class data frame containing the output of the differential splicing analysis.
#'
#' @importFrom plyr join
#' @importFrom stats ks.test na.omit p.adjust p.adjust.methods t.test wilcox.test
#' @import methods
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' # Define cell groups for analysis
#' df.pheno <- marvel.demo$SplicePheno
#' cell.group.g1 <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
#' cell.group.g2 <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]
#'
#' # DE
#' results <- CompareValues.PSI(MarvelObject=marvel.demo,
#'                              cell.group.g1=cell.group.g1,
#'                              cell.group.g2=cell.group.g2,
#'                              min.cells=5,
#'                              method="t.test",
#'                              method.adjust="fdr",
#'                              event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
#'                              show.progress=FALSE
#'                              )
#'
#' # Check output
#' head(results)

CompareValues.PSI <- function(MarvelObject, cell.group.g1, cell.group.g2, downsample=FALSE,  seed=1, min.cells=25, pct.cells=NULL, method, nboots=1000, n.permutations=1000, method.adjust="fdr", event.type, show.progress=TRUE, annotate.outliers=TRUE, n.cells.outliers=10, assign.modality=TRUE) {

    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.group.g1 <- cell.group.g1
    cell.group.g2 <- cell.group.g2
    downsample <- downsample
    min.cells <- min.cells
    pct.cells <- pct.cells
    method <- method
    nboots <- nboots
    n.permutations <- n.permutations
    method.adjust <- method.adjust
    event.type <- event.type
    show.progress <- show.progress
    annotate.outliers <- annotate.outliers
    n.cells.outliers <- n.cells.outliers
    assign.modality <- assign.modality
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- do.call(rbind.data.frame, MarvelObject$PSI)
    #df.pheno <- MarvelObject$SplicePheno
    #df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    #cell.group.g1 <- cell.group.g1
    #cell.group.g2 <- cell.group.g2
    #downsample <- FALSE
    #min.cells <- 10
    #pct.cells <- NULL
    #method <- "ad"
    #method.adjust <- "fdr"
    #event.type <- c("A3SS")
    #show.progress <- TRUE
    #annotate.outliers <- TRUE
    #n.cells.outliers <- 5
    #assign.modality <- TRUE
    
    ############################################################
    
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Subset relevant events
    df.feature <- df.feature[which(df.feature$event_type %in% event.type), ]
    df <- df[df.feature$tran_id,]
    
    # Retrieve sample ids
        # Group 1
        sample.ids.1 <- cell.group.g1
        
        # Group 2
        sample.ids.2 <- cell.group.g2
        
        # Subset
        df <- df[, c(sample.ids.1, sample.ids.2)]
    
    # Downsample
    if(downsample==TRUE) {
        
        # Retrieve lowest denominator
        n.cells.downsample <- min(length(sample.ids.1), length(sample.ids.2))
        
        # Downsample
        set.seed(seed)
        sample.ids.1 <- sample(sample.ids.1, size=n.cells.downsample, replace=FALSE)
        sample.ids.2 <- sample(sample.ids.2, size=n.cells.downsample, replace=FALSE)
        
        # Subset cells
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% c(sample.ids.1, sample.ids.2)), ]
        df <- df[, df.pheno$sample.id]
        
        # Track progress
        message(paste(length(sample.ids.1), " cells found in Group 1", sep=""))
        message(paste(length(sample.ids.2), " cells found in Group 2", sep=""))
        message(paste("Both Group 1 and 2 downsampled to ", n.cells.downsample, " cells", sep=""))
        
    }
    
    # Subset events with sufficient cells
        # Group 1
        df.small <- df[, which(names(df) %in% sample.ids.1)]
        . <- apply(df.small, 1, function(x) {sum(!is.na(x))})
        
        if(is.null(pct.cells)) {
            
            tran_ids.1 <- names(.)[which(. >= min.cells)]
        
        } else {
            
            . <- ./nrow(df.pheno) * 100
            
            tran_ids.1 <- names(.)[which(. >= pct.cells)]
            
        }
        
        # Group 2
        df.small <- df[, which(names(df) %in% sample.ids.2)]
        . <- apply(df.small, 1, function(x) {sum(!is.na(x))})
        
        
        if(is.null(pct.cells)) {
            
            tran_ids.2 <- names(.)[which(. >= min.cells)]
        
        } else {
            
            . <- ./nrow(df.pheno) * 100
            
            tran_ids.2 <- names(.)[which(. >= pct.cells)]
            
        }
       
        # Subset overlaps
        overlap <- intersect(tran_ids.1, tran_ids.2)
        df.feature <- df.feature[which(df.feature$tran_id %in% overlap), ]
        df <- df[overlap, ]
        
        # Report progress
        message(paste(length(tran_ids.1), " expressed events identified in Group 1", sep=""))
        message(paste(length(tran_ids.2), " expressed events identified in Group 2", sep=""))
        message(paste(length(overlap), " expressed events identified in BOTH Group 1 and Group 2", sep=""))
        
    # Remove events in which PSI values are constant across all cells
    if(method=="permutation") {
        
        # Retrieve non-variable events
        . <- apply(df, 1, function(x) {length(unique(x[!is.na(x)]))})
        index.rm <- which(.==1)
        
        # Remove non-variable events
        if(length(index.rm) != 0) {
            
            df <- df[-index.rm, ]
            
        }
        
        message(paste(length(index.rm), " non-variable events identified and removed", sep=""))
        message(paste(nrow(df), " events retained", sep=""))
        
    }
    
    ######################################################################
    
    # Statistical test
    tran_ids <- row.names(df)
    
    n.cells.x <- NULL
    n.cells.y <- NULL
    mean.x <- NULL
    mean.y <- NULL
    mean.diff <- NULL
    p.val <- NULL
    statistic <- NULL
    
    if(show.progress==TRUE) {
        
        pb <- txtProgressBar(1, length(tran_ids), style=3)
        
    }

    if(method != "permutation") {
        
    ######################################################################
    ############### NON-PERMUTATION APPROACH (OTHER THAN DTS) ############
    ######################################################################
    
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
                
            } else if(method=="ad") {
                
                error.check <- tryCatch(kSamples::ad.test(x, y), error=function(err) "Error")
                
                if(error.check[1] == "Error") {
                    
                    statistic[i] <- 0
                    p.val[i] <- 1
                    
                } else {
                    
                    stats <- kSamples::ad.test(x, y, method="asymptotic")$ad
                    
                    statistic[i] <- stats[1,1]
                    p.val[i] <- stats[1,3]
                    
                }
                
            } else if(method=="dts"){
                
                stats <- twosamples::dts_test(x, y, nboots=nboots)
                
                statistic[i] <- stats[1]
                p.val[i] <- stats[2]
                
                
            }
                    
            # Track progress
            if(show.progress==TRUE) {
           
                #setTxtProgressBar(pb, i)
                setTxtProgressBar(pb, i)
                #message(i)
                
            }
            
        }
        
    ######################################################################
    ######################### PERMUTATION APPROACH #######################
    ######################################################################
    
    } else if(method=="permutation") {
    
        # Compute p-values
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
            
            # Compute p-value
                # Compute observed mean
                mean.diff.obs <-  mean(y) - mean(x)
                
                # Subset event
                df.small <- df[tran_ids[i], ]
                df.small <- na.omit(df.small)
                            
                # Build null distribution
                set.seed(seed)
                
                mean.diff.perm <- NULL
                
                for(j in 1:n.permutations) {
                    
                    # Shuffle group labels
                    names(df.small) <- sample(x=names(df.small),
                                              size=length(names(df.small)),
                                              replace=FALSE,
                                              )
                    
                    # Retrieve values
                    x.perm <- as.numeric(df.small[, intersect(names(df.small), sample.ids.1)])
                    y.perm <- as.numeric(df.small[, intersect(names(df.small), sample.ids.2)])
                    
                    # Compute mean
                    mean.diff.perm[j] <- mean(y.perm) - mean(x.perm)
                                    
                }
                
                # Compute p-value
                if(mean.diff.obs < 0) {
                    
                    statistic[i] <- NA
                    p.val[i] <- sum(mean.diff.perm < mean.diff.obs)/length(mean.diff.perm)
                    
                } else if(mean.diff.obs >= 0){
                    
                    statistic[i] <- NA
                    p.val[i] <- sum(mean.diff.perm > mean.diff.obs)/length(mean.diff.perm)
                    
                }
                
            # Track progress
            if(show.progress==TRUE) {
           
                #setTxtProgressBar(pb, i)
                setTxtProgressBar(pb, i)
                #message(i)
                
            }

        }
        
    }
    
    ######################################################################
    
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
    #results <- results[which(!is.na(results$p.val)), ]
    results$p.val[which(is.na(results$p.val))] <- 1
    #results$statistic[which(results$p.val < 0)]  <- results$statistic[which(results$p.val < 0)] * -1
    #results$p.val <- abs(results$p.val)
    results <- results[order(results$p.val), ]
    #results <- results[order(results$statistic, decreasing=TRUE), ]
    
    # Adjust for multiple testing
    if(method != "dts" & method != "permutation") {
        
        results$p.val.adj <- p.adjust(results$p.val, method=method.adjust, n=length(results$p.val))
        #results$p.val <- NULL
        
    } else if(method=="dts" | method == "permutation"){
    
        results$p.val.adj <- results$p.val
        #results$p.val <- NULL
        
    }
    
    # Annotate with feature metadata
    results <- join(results, df.feature, by="tran_id", type="left")
    cols.1 <- names(df.feature)
    cols.2 <- setdiff(names(results), names(df.feature))
    results <- results[,c(cols.1, cols.2)]
    
    #############################################################################
    ########################### ASSIGN MODALITIES  ##############################
    #############################################################################
    
    if(assign.modality==TRUE) {
        
        # Track progress
        message("Assigning modalities...")
        
        # Assign modality
            # Group 1
            modality <- AssignModality(MarvelObject=MarvelObject,
                                       sample.ids=cell.group.g1,
                                       min.cells=min.cells,
                                       bimodal.adjust=TRUE,
                                       seed=seed,
                                       tran_ids=results$tran_id
                                       )
                                
           modality.g1 <- modality$Modality$Results
           
           # Group 2
           modality <- AssignModality(MarvelObject=MarvelObject,
                                      sample.ids=cell.group.g2,
                                      min.cells=min.cells,
                                      bimodal.adjust=TRUE,
                                      seed=seed,
                                      tran_ids=results$tran_id
                                      )
                               
           modality.g2 <- modality$Modality$Results
           
        # Annotate result table
            # Annotate
            names(modality.g1)[which(names(modality.g1)=="modality.bimodal.adj")] <- "modality.bimodal.adj.g1"
            names(modality.g2)[which(names(modality.g2)=="modality.bimodal.adj")] <- "modality.bimodal.adj.g2"
            results <- join(results, modality.g1[,c("tran_id", "modality.bimodal.adj.g1")], by="tran_id", type="left")
            results <- join(results, modality.g2[,c("tran_id", "modality.bimodal.adj.g2")], by="tran_id", type="left")

            # Check for non-matches
            sum(is.na(results$modality.bimodal.adj.g1))
            sum(is.na(results$modality.bimodal.adj.g2))

    }
    
    #############################################################################
    ########################## ANNOTATE FALSE +VE'S  ############################
    #############################################################################
    
    if(annotate.outliers==TRUE) {
        
        # Track progress
        message("Identifying outliers...")
        
        # Check num of cells !=1 for included -> included modality changes
            # Subset relevant modality changes
            index.1 <- grep("Included", results$modality.bimodal.adj.g1, fixed=TRUE)
            index.2 <- grep("Included", results$modality.bimodal.adj.g2, fixed=TRUE)
            index <- intersect(index.1, index.2)
            results.small <- results[index, ]
            
            if(length(index) != 0) {
            
                # Group 1
                    # Subset matrix
                    df.small <- df[results.small$tran_id, sample.ids.1]
                    
                    # Compute n cells
                    . <- apply(df.small, 1, function(x) {sum(x!=1, na.rm=TRUE)})
                    . <- data.frame("tran_id"=names(.), "n.cells.outliers.g1"=as.numeric(.), stringsAsFactors=FALSE)
                    
                    # Annotate result table
                    results.small <- join(results.small, ., by="tran_id", type="left")
                    
                # Group 2
                    # Subset matrix
                    df.small <- df[results.small$tran_id, sample.ids.2]
                    
                    # Compute n cells
                    . <- apply(df.small, 1, function(x) {sum(x!=1, na.rm=TRUE)})
                    . <- data.frame("tran_id"=names(.), "n.cells.outliers.g2"=as.numeric(.), stringsAsFactors=FALSE)
                    
                    # Annotate result table
                    results.small <- join(results.small, ., by="tran_id", type="left")
                    
                # Re-code p-values below threshold
                index <- which(results.small$n.cells.outliers.g1 < n.cells.outliers & results.small$n.cells.outliers.g2 < n.cells.outliers)
                
                #results.small$p.val.adj[index] <- 1
                
                # Indicate if outliers were adjusted for
                results.small$outliers <- FALSE
                results.small$outliers[index] <- TRUE
                
                # Save as new object
                results.small.included <- results.small
                
            } else {
                
                results.small.included <- NULL
                
            }
        
        # Check num of cells !=0 for excluded -> excluded modality changes
            # Subset relevant modality changes
            index.1 <- grep("Excluded", results$modality.bimodal.adj.g1, fixed=TRUE)
            index.2 <- grep("Excluded", results$modality.bimodal.adj.g2, fixed=TRUE)
            index <- intersect(index.1, index.2)
            results.small <- results[index, ]
            
            if(length(index) != 0) {
            
                # Group 1
                    # Subset matrix
                    df.small <- df[results.small$tran_id, sample.ids.1]
                    
                    # Compute n cells
                    . <- apply(df.small, 1, function(x) {sum(x!=0, na.rm=TRUE)})
                    . <- data.frame("tran_id"=names(.), "n.cells.outliers.g1"=as.numeric(.), stringsAsFactors=FALSE)
                    
                    # Annotate result table
                    results.small <- join(results.small, ., by="tran_id", type="left")
                    
                # Group 2
                    # Subset matrix
                    df.small <- df[results.small$tran_id, sample.ids.2]
                    
                    # Compute n cells
                    . <- apply(df.small, 1, function(x) {sum(x!=0, na.rm=TRUE)})
                    . <- data.frame("tran_id"=names(.), "n.cells.outliers.g2"=as.numeric(.), stringsAsFactors=FALSE)
                    
                    # Annotate result table
                    results.small <- join(results.small, ., by="tran_id", type="left")
                    
                # Re-code p-values below threshold
                index <- which(results.small$n.cells.outliers.g1 < n.cells.outliers & results.small$n.cells.outliers.g2 < n.cells.outliers)
                
                #results.small$p.val.adj[index] <- 1
                
                # Indicate if outliers
                results.small$outliers <- FALSE
                results.small$outliers[index] <- TRUE
                
                # Save as new object
                results.small.excluded <- results.small
                
            } else {
                
                results.small.excluded <- NULL
                
            }
            
        # Merge
        if(!is.null(results.small.included) | !is.null(results.small.excluded)) {
            
            if(is.null(results.small.included) & nrow(results.small.excluded) == nrow(results)) {
                    
                results <- results.small.excluded
                
            } else if(is.null(results.small.included) & nrow(results.small.excluded) != nrow(results)) {
                
                tran_ids <- setdiff(results$tran_id, results.small.excluded$tran_id)
                results.small.. <- results[which(results$tran_id %in% tran_ids), ]
                results.small..$n.cells.outliers.g1 <- 0
                results.small..$n.cells.outliers.g2 <- 0
                results.small..$outliers <- FALSE
                
                results <- rbind.data.frame(results.small.., results.small.excluded)
            
            } else if(nrow(results.small.included) == nrow(results) & is.null(results.small.excluded)) {
                
                results <- results.small.included
                
            } else if(nrow(results.small.included) == nrow(results) & is.null(results.small.excluded)) {
                
                tran_ids <- setdiff(results$tran_id, results.small.included$tran_id)
                results.small.. <- results[which(results$tran_id %in% tran_ids), ]
                results.small..$n.cells.outliers.g1 <- 0
                results.small..$n.cells.outliers.g2 <- 0
                results.small..$outliers <- FALSE
                
                results <- rbind.data.frame(results.small.., results.small.included)
                
            } else if(!is.null(results.small.included) & !is.null(results.small.excluded)){
                
                # Retrieve events not needing outlier adjustment
                tran_ids <- setdiff(results$tran_id, c(results.small.included$tran_id, results.small.excluded$tran_id))
            
                results.small.. <- results[which(results$tran_id %in% tran_ids), ]
                
                # Match columns
                results.small..$n.cells.outliers.g1 <- 0
                results.small..$n.cells.outliers.g2 <- 0
                results.small..$outliers <- FALSE
                
                # Merge
                results <- rbind.data.frame(results.small.., results.small.included, results.small.excluded)
            
                # Reorder by p-values
                results <- results[order(results$p.val), ]
                
            }
            
        } else {
            
            # Match columns
            results$n.cells.outliers.g1 <- 0
            results$n.cells.outliers.g2 <- 0
            results$outliers <- FALSE
            
        }
        
    } else {
        
            # Match columns
            results$n.cells.outliers.g1 <- 0
            results$n.cells.outliers.g2 <- 0
            results$outliers <- FALSE
        
    }
        
    
    #############################################################################
    
    # Convert fraction to %
    results$mean.g1 <- results$mean.g1 * 100
    results$mean.g2 <- results$mean.g2 * 100
    results$mean.diff <- results$mean.diff * 100
    
    # Report result summary
    #message(paste(sum(results$p.val.adj < 0.10), " DE splicing events < 0.10 adjusted p-value", sep=""))
    #message(paste(sum(results$p.val.adj < 0.05), " DE splicing events < 0.05 adjusted p-value", sep=""))
    #message(paste(sum(results$p.val.adj < 0.01), " DE splicing events < 0.01 adjusted p-value", sep=""))
 
    # Return result table (not new MARVEL object)
    return(results)
        
}
