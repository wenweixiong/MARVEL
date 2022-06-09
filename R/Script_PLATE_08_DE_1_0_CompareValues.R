#' @title Differential splicing and gene expression analysis
#'
#' @description Performs differential splicing and gene expression analysis between 2 groups of cells. This is a wrapper function for \code{CompareValues.PSI} and \code{CompareValues.Exp} functions.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2.
#' @param downsample Logical value. If set to \code{TRUE}, the number of cells in each cell group will be downsampled to the sample size of the smaller cell group so that both cell groups will have the sample size prior to differential expression analysis. Default is \code{FALSE}.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event or genes for the event or genes to be included for differential splicing analysis.
#' @param pct.cells Numeric value. The minimum percentage of cells expressing the splicing event or genes for the event or genes to be included for differential splicing analysis. If \code{pct.cells} is specified, then \code{pct.cells} will be used as threshold instead of \code{min.cells}.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"kuiper"}, \code{"ad"}, \code{"dts"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Kuiper, Anderson-Darling, DTS, Wilcox, and t-test, respectively.
#' @param nboots Numeric value.Only applicable when \code{level} set to \code{"splicing"}.  When \code{method} set to \code{"dts"}, the number of bootstrap iterations for computing the p-value.
#' @param n.permutations Numeric value. Only applicable when \code{level} set to \code{"splicing"}. When \code{method} set to \code{"permutation"}, this argument indicates the number of permutations to perform for generating the null distribution for subsequent p-value inference. Default is \code{1000} times.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for differential splicing or gene expression analysis, respectively.
#' @param event.type Character string. Only applicable when \code{level} set to \code{"splicing"}. Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param show.progress Logical value. If set to \code{TRUE}, progress bar will be displayed so that users can estimate the time needed for differential analysis. Default value is \code{TRUE}.
#' @param annotate.outliers Numeric value. Only applicable when \code{level} set to \code{"splicing"}. When set to \code{TRUE}, statistical difference in PSI values between the two cell groups that is driven by outlier cells will be annotated.
#' @param n.cells.outliers Numeric value. Only applicable when \code{level} set to \code{"splicing"}. When \code{method} set to \code{"dts"}, the minimum number of cells with non-1 or non-0 PSI values for included-to-included or excluded-to-excluded modality change, respectively. The p-values will be re-coded to 1 when both cell groups have less than this minimum number of cells. This is to avoid false positive results.
#' @param assign.modality Logical value. Only applicable when \code{level} set to \code{"splicing"}. If set to \code{TRUE} (default), modalities will be assigned to each cell group.
#'
#' @return An object of class S3 containing with new slot \code{MarvelObject$DE$PSI$Table[["method"]]} or \code{MarvelObject$DE$Exp$Table} when \code{level} option specified as \code{"splicing"} or \code{"gene"}, respectively.
#'
#' @importFrom plyr join
#' @importFrom waddR wasserstein.test
#' @importFrom kSamples ad.test
#' @import stats
#' @import methods
#' @import kuiper.2samp
#' @import twosamples
#' @import scran
#'
#' @export

CompareValues <- function(MarvelObject, cell.group.g1, cell.group.g2, downsample=FALSE, min.cells=25, pct.cells=NULL, method, nboots=1000, n.permutations=1000, method.adjust="fdr", level, event.type=NULL, show.progress=TRUE, annotate.outliers=TRUE, n.cells.outliers=10, assign.modality=TRUE) {

    # Define arguments
    MarvelObject <- MarvelObject
    cell.group.g1 <- cell.group.g1
    cell.group.g2 <- cell.group.g2
    downsample <- downsample
    min.cells <- min.cells
    pct.cells <- pct.cells
    method <- method
    nboots <- nboots
    n.permutations <- n.permutations
    method.adjust <- method.adjust
    level <- level
    event.type <- event.type
    show.progress <- show.progress
    annotate.outliers <- annotate.outliers
    n.cells.outliers <- n.cells.outliers
    assign.modality <- assign.modality
    
    # Example arguments (splicing)
    #MarvelObject <- marvel
    #df <- do.call(rbind.data.frame, MarvelObject$PSI)
    #df.pheno <- MarvelObject$SplicePheno
    #df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    #cell.group.g1 <- cell.group.g1
    #cell.group.g2 <- cell.group.g2
    #downsample <- FALSE
    #min.cells <- 25
    #pct.cells <- NULL
    #method <- c("ad", "dts")
    #method.adjust <- "fdr"
    #event.type <- c("ALE")
    #show.progress <- TRUE
    #nboots <- 2000
    #annotate.outliers <- TRUE
    #n.cells.outliers <- 10
    
    if(level=="splicing") {
        
        for(i in 1:length(method)) {
            
            # Track progress
            print(paste("Running ", method[i], " ...", sep=""))
            
            results <- CompareValues.PSI(MarvelObject=MarvelObject,
                                         cell.group.g1=cell.group.g1,
                                         cell.group.g2=cell.group.g2,
                                         downsample=downsample,
                                         min.cells=min.cells,
                                         method=method[i],
                                         nboots=nboots,
                                         n.permutations=n.permutations,
                                         method.adjust=method.adjust,
                                         event.type=event.type,
                                         annotate.outliers=annotate.outliers,
                                         n.cells.outliers=n.cells.outliers,
                                         show.progress=show.progress,
                                         assign.modality=assign.modality
                                         )
                        
            # Save to new slot
            MarvelObject$DE$PSI$Table[[method[i]]] <- results
            
        }
        
        # Return final object
        return(MarvelObject)
                                     
    } else if(level=="gene") {
        
        CompareValues.Exp(MarvelObject=MarvelObject,
                          cell.group.g1=cell.group.g1,
                          cell.group.g2=cell.group.g2,
                          downsample=downsample,
                          min.cells=min.cells,
                          pct.cells=pct.cells,
                          method=method,
                          method.adjust=method.adjust,
                          nboots=nboots,
                          show.progress=show.progress
                          )
        
    }
    
}
