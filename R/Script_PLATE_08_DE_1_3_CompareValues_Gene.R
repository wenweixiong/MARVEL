#' @title Differential gene expression analysis
#'
#' @description Performs differential gene expression analysis between 2 groups of cells.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2.
#' @param downsample Logical value. If set to \code{TRUE}, the number of cells in each cell group will be downsampled to the sample size of the smaller cell group so that both cell groups will have the sample size prior to differential expression analysis. Default is \code{FALSE}.
#' @param seed Numeric value. The seed number for the random number generator to ensure reproducibility during during down-sampling of cells when \code{downsample} set to \code{TRUE}.
#' @param min.cells Numeric value. The minimum no. of cells expressing the gene for the gene to be included for differential splicing analysis.
#' @param pct.cells Numeric value. The minimum no. of cells expressing the gene for the gene to be included for differential splicing analysis. If \code{pct.cells} is specified, then \code{pct.cells} will be used as threshold instead of \code{min.cells}.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"kuiper"}, \code{"ad"}, \code{"dts"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Kuiper, Anderson-Darling, DTS, Wilcox, and t-test, respectively. Additional option is \code{"mast"}. If set to \code{"mast"} is specified, the log2fc and p-values will be corrected using the gene detection rate as per the \code{MAST} package tutorial.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param show.progress Logical value. If set to \code{TRUE}, progress bar will be displayed so that users can estimate the time needed for differential analysis. Default value is \code{TRUE}.
#' @param nboots Numeric value. When \code{method} set to \code{"dts"}, the number of bootstrap iterations for computing the p-value.
#' @param custom.gene_ids Character string. Instead of specified the genes to include for DE analysis with \code{min.cells}, users may input a custom vector of gene IDs to include for DE analysis.
#' @param mast.method Character string. As per the \code{method} option of the \code{zlm} function from the \code{MAST} package. Default is \code{"bayesglm"}, other options are \code{"glm"} and \code{"glmer"}.
#' @param mast.ebayes Logical value. As per the \code{ebayes} option of the \code{zlm} function from the \code{MAST} package. Default is \code{TRUE}.
#'
#' @return An object of class S3 new slot \code{MarvelObject$DE$Exp$Table}.
#'
#' @importFrom plyr join
#' @importFrom stats ks.test na.omit p.adjust p.adjust.methods t.test wilcox.test
#' @import methods
#' @importFrom utils txtProgressBar setTxtProgressBar
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
#' marvel.demo <- CompareValues.Exp(MarvelObject=marvel.demo,
#'                                  cell.group.g1=cell.group.g1,
#'                                  cell.group.g2=cell.group.g2,
#'                                  min.cells=5,
#'                                  method="t.test",
#'                                  method.adjust="fdr",
#'                                  show.progress=FALSE
#'                                  )
#'
#' # Check output
#' head(marvel.demo$DE$Exp$Table)

CompareValues.Exp <- function(MarvelObject, cell.group.g1=NULL, cell.group.g2=NULL, downsample=FALSE, seed=1, min.cells=25, pct.cells=NULL, method, method.adjust, show.progress=TRUE, nboots=1000, custom.gene_ids=NULL, mast.method="bayesglm", mast.ebayes=TRUE) {

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
    
    # Define arguments
    #custom.gene_ids <- df$gene_id
    #df <- marvel$Exp
    #df.pheno <- marvel$SplicePheno
    #df.feature <- marvel$GeneFeature
    #cell.group.g1 <- cell.group.g1
    #cell.group.g2 <- cell.group.g2
    #downsample <- FALSE
    #min.cells <- 3
    #pct.cells <- NULL
    #method <- "mast"
    #method.adjust <- "fdr"
    #show.progress <- TRUE
    #mast.method <- "bayesglm"
    #mast.ebayes <- TRUE
    
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
  
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
        set.seed(seed)
        sample.ids.1.small <- sample(sample.ids.1, size=n.cells.downsample, replace=FALSE)
        sample.ids.2.small <- sample(sample.ids.2, size=n.cells.downsample, replace=FALSE)
        
        # Subset cells
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% c(sample.ids.1.small, sample.ids.2.small)), ]
        df <- df[, df.pheno$sample.id]
        
        # Track progress
        message(paste(length(sample.ids.1), " cells found in Group 1", sep=""))
        message(paste(length(sample.ids.2), " cells found in Group 2", sep=""))
        message(paste("Both Group 1 and 2 downsampled to ", n.cells.downsample, " cells", sep=""))
        
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
        message(paste(length(gene_ids.1), " expressed genes identified in Group 1", sep=""))
        message(paste(length(gene_ids.2), " expressed genes identified in Group 2", sep=""))
        message(paste(length(all), " expressed genes identified in EITHER Group 1 or Group 2", sep=""))
    
    } else {
        
        df.feature <- df.feature[which(df.feature$gene_id %in% custom.gene_ids), ]
        df <- df[df.feature$gene_id, ]
        message(paste(length(custom.gene_ids), " custom genes specified", sep=""))
        
    }
    
    ########################################################
    
    # Statistical test
    if(method != "mast") {
        
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
                
            } else if(method=="ad") {
                
                error.check <- tryCatch(kSamples::ad.test(x, y), error=function(err) "Error")
                
                if(error.check[1] == "Error") {
                    
                    statistic[i] <- 0
                    p.val[i] <- 1
                    
                } else {
                    
                    . <- kSamples::ad.test(x, y, method="asymptotic")$ad
                    
                    statistic[i] <- .[1,1]
                    p.val[i] <- .[1,3]
                    
                }
                
            } else if(method=="dts"){
                
                . <- twosamples::dts_test(x, y, nboots=nboots)
                
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
                              
    } else if(method=="mast"){
        
        # Prepare cdata
            # Indicate group 1,2
            cdata <- data.frame("sample.id"=c(sample.ids.1, sample.ids.2))
            cdata$condition <- ifelse(cdata$sample.id %in% sample.ids.1, "g1", "g2")
            cdata$condition <- factor(cdata$condition, levels=c("g1", "g2"))
            
            # Retrieve original GE matrix
            df.exp.master <- MarvelObject$Exp
            row.names(df.exp.master) <- df.exp.master$gene_id
            df.exp.master$gene_id <- NULL
            df.exp.master <- df.exp.master[,cdata$sample.id]
            
            # Compute n expressed genes
            . <- apply(df.exp.master, 2, function(x) {sum(x != 0)})
            . <- data.frame("sample.id"=names(.),
                            "n.genes"=as.numeric(.),
                            stringsAsFactors=FALSE
                            )
                            
            # Compute gene detection rate
            .$cngeneson <- scale(.$n.genes)
            cdata <- join(cdata, ., by="sample.id", type="left")
            
            # Format for MAST
            names(cdata)[which(names(cdata)=="sample.id")] <- "wellKey"
            row.names(cdata) <- cdata$wellKey
            
        # Prepare fdata
        fdata <- df.feature
        names(fdata)[which(names(fdata)=="gene_id")] <- "primerid"
        row.names(fdata) <- fdata$primerid
        
        # Prepare SingleCellAssay object
        df <- df[, cdata$wellKey]
        df <- df[fdata$primerid, ]
        sca <- MAST::FromMatrix(as.matrix(df), cdata, fdata)
        
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
        
        # Match MARVEL format
            # Subset relevant columns
            results <- fcHurdle
            results <- results[,c("primerid", "Pr(>Chisq)", "coef")]
            names(results) <- c("gene_id", "p.val", "log2fc")
            
            # Compute other metrices
            gene_ids <- df.feature$gene_id
            
            n.cells.x <- NULL
            n.cells.y <- NULL
            mean.x <- NULL
            mean.y <- NULL

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
                
                # Track progress
                if(show.progress==TRUE) {
                    
                    setTxtProgressBar(pb, i)
                    
                }
                
            }
            
            results.2 <- data.frame("gene_id"=gene_ids,
                                    "n.cells.g1"=n.cells.x, "n.cells.g2"=n.cells.y,
                                    "mean.g1"=mean.x, "mean.g2"=mean.y,
                                    "statistic"=NA,
                                    stringsAsFactors=FALSE
                                    )
            
            results <- join(results, results.2, by="gene_id", type="left")
            
            # Match column order
            cols <- c("gene_id",
                      "n.cells.g1", "n.cells.g2", "mean.g1", "mean.g2",
                      "log2fc", "statistic", "p.val"
                      )
            results <- results[,cols]
        
    }
    
    ########################################################

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
    #message(paste(sum(results$p.val.adj < 0.10), " DE genes < 0.10 adjusted p-value", sep=""))
    #message(paste(sum(results$p.val.adj < 0.05), " DE genes < 0.05 adjusted p-value", sep=""))
    #message(paste(sum(results$p.val.adj < 0.01), " DE genes < 0.01 adjusted p-value", sep=""))
    
    # Save to new slot
    if(is.null(custom.gene_ids[1])) {
        
        MarvelObject$DE$Exp$Table <- results
        
    } else {
        
        MarvelObject$DE$Exp.Custom$Table <- results
        
    }
  
    return(MarvelObject)
        
}
