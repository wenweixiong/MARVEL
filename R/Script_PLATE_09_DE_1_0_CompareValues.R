#' @title Differential splicing and gene expression analysis
#'
#' @description Performs differential splicing and gene expression analysis between 2 groups of cells. This is a wrapper function for \code{CompareValues.PSI} and \code{CompareValues.Exp} functions.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2.
#' @param downsample Logical value. If set to \code{TRUE}, the number of cells in each cell group will be downsampled to the sample size of the smaller cell group so that both cell groups will have the sample size prior to differential expression analysis. Default is \code{FALSE}.
#' @param seed Numeric value. The seed number for the random number generator to ensure reproducibility during during down-sampling of cells when \code{downsample} set to \code{TRUE}.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event or genes for the event or genes to be included for differential splicing analysis.
#' @param pct.cells Numeric value. The minimum percentage of cells expressing the splicing event or genes for the event or genes to be included for differential splicing analysis. If \code{pct.cells} is specified, then \code{pct.cells} will be used as threshold instead of \code{min.cells}.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"kuiper"}, \code{"ad"}, \code{"dts"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Kuiper, Anderson-Darling, DTS, Wilcox, and t-test, respectively. Additional \code{"mast"} option is available for differential gene expression analysis. If \code{"mast"} is specified, the log2fc and p-values will be corrected using the gene detection rate as per the \code{MAST} package tutorial.
#' @param nboots Numeric value.Only applicable when \code{level} set to \code{"splicing"}.  When \code{method} set to \code{"dts"}, the number of bootstrap iterations for computing the p-value.
#' @param n.permutations Numeric value. Only applicable when \code{level} set to \code{"splicing"}. When \code{method} set to \code{"permutation"}, this argument indicates the number of permutations to perform for generating the null distribution for subsequent p-value inference. Default is \code{1000} times.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for differential splicing or gene expression analysis, respectively.
#' @param event.type Character string. Only applicable when \code{level} set to \code{"splicing"}. Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param show.progress Logical value. If set to \code{TRUE}, progress bar will be displayed so that users can estimate the time needed for differential analysis. Default value is \code{TRUE}.
#' @param annotate.outliers Numeric value. Only applicable when \code{level} set to \code{"splicing"}. When set to \code{TRUE}, statistical difference in PSI values between the two cell groups that is driven by outlier cells will be annotated.
#' @param n.cells.outliers Numeric value. Only applicable when \code{level} set to \code{"splicing"}. When \code{method} set to \code{"dts"}, the minimum number of cells with non-1 or non-0 PSI values for included-to-included or excluded-to-excluded modality change, respectively. The p-values will be re-coded to 1 when both cell groups have less than this minimum number of cells. This is to avoid false positive results.
#' @param assign.modality Logical value. Only applicable when \code{level} set to \code{"splicing"}. If set to \code{TRUE} (default), modalities will be assigned to each cell group.
#' @param custom.gene_ids Character string. Only applicable when \code{level} set to \code{"gene"}. Instead of specified the genes to include for DE analysis with \code{min.cells}, users may input a custom vector of gene IDs to include for DE analysis.
#' @param psi.method Vector of character string(s). Only applicable when \code{level} set to \code{"gene.spliced"} and when \code{CompareValues} function has been ran with \code{level} set to \code{"splicing"} earlier. To include significant events from these method(s) for differential gene expression analysis.
#' @param psi.pval Vector of numeric value(s). Only applicable when \code{level} set to \code{"gene.spliced"} and when \code{CompareValues} function has been ran with \code{level} set to \code{"splicing"} earlier. The adjusted p-value, below which, the splicing event is considered differentially spliced, and the corresponding genes will be included for differential gene expression analysis.
#' @param psi.delta Numeric value. Only applicable when \code{level} set to \code{"gene.spliced"} and when \code{CompareValues} function has been ran with \code{level} set to \code{"splicing"} earlier. The absolute difference in mean PSI values between \code{cell.group.g1} and \code{cell.group.g1}, above which, the splicing event is considered differentially spliced, and the corresponding genes will be included for differential gene expression analysis.
#' @param method.de.gene Character string. Only applicable when \code{level} set to \code{"gene.spliced"} and when \code{CompareValues} function has been ran with \code{level} set to \code{"splicing"} earlier. Same as \code{method}.
#' @param method.adjust.de.gene Character string. Only applicable when \code{level} set to \code{"gene.spliced"} and when \code{CompareValues} function has been ran with \code{level} set to \code{"splicing"} earlier. Same as \code{method.adjust}.
#' @param mast.method Character string. Only applicable when \code{level} set to \code{"gene"} or \code{"gene.spliced"}. As per the \code{method} option of the \code{zlm} function from the \code{MAST} package. Default is \code{"bayesglm"}, other options are \code{"glm"} and \code{"glmer"}.
#' @param mast.ebayes Logical value. Only applicable when \code{level} set to \code{"gene"} or \code{"gene.spliced"}. As per the \code{ebayes} option of the \code{zlm} function from the \code{MAST} package. Default is \code{TRUE}.
#' @param seed Numeric value. The seed number for the random number generator to ensure reproducibility during during down-sampling of cells when \code{downsample} set to \code{TRUE}, during permutation testing when \code{method} set to \code{"permutation"}, and during modality assignment which will be performed automatically.
#'
#' @return An object of class S3 containing with new slot \code{MarvelObject$DE$PSI$Table[["method"]]} or \code{MarvelObject$DE$Exp$Table} when \code{level} option specified as \code{"splicing"} or \code{"gene"}, respectively.
#'
#' @importFrom plyr join
#' @import methods
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
#' marvel.demo <- CompareValues(MarvelObject=marvel.demo,
#'                              cell.group.g1=cell.group.g1,
#'                              cell.group.g2=cell.group.g2,
#'                              min.cells=5,
#'                              method="t.test",
#'                              method.adjust="fdr",
#'                              level="splicing",
#'                              event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
#'                              show.progress=FALSE
#'                              )
#'
#' # Check output
#' head(marvel.demo$DE$PSI$Table[["ad"]])

CompareValues <- function(MarvelObject, cell.group.g1=NULL, cell.group.g2=NULL, downsample=FALSE, seed=1, min.cells=25, pct.cells=NULL, method=NULL, nboots=1000, n.permutations=1000, method.adjust="fdr", level, event.type=NULL, show.progress=TRUE, annotate.outliers=TRUE, n.cells.outliers=10, assign.modality=TRUE, custom.gene_ids=NULL, psi.method=NULL, psi.pval=NULL, psi.delta=NULL, method.de.gene=NULL, method.adjust.de.gene=NULL, mast.method="bayesglm", mast.ebayes=TRUE) {

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
    custom.gene_ids <- custom.gene_ids
    psi.method <- psi.method
    psi.pval <- psi.pval
    psi.delta <- psi.delta
    method.de.gene <- method.de.gene
    method.adjust.de.gene <-  method.adjust.de.gene
    
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
            message(paste("Running ", method[i], " ...", sep=""))
            
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
                                         assign.modality=assign.modality,
                                         seed=seed
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
                          seed=seed,
                          min.cells=min.cells,
                          pct.cells=pct.cells,
                          method=method,
                          method.adjust=method.adjust,
                          nboots=nboots,
                          show.progress=show.progress,
                          custom.gene_ids=custom.gene_ids,
                          mast.method=mast.method,
                          mast.ebayes=mast.ebayes
                          )
        
    } else if(level=="gene.spliced") {
        
        CompareValues.Exp.Spliced(MarvelObject=MarvelObject,
                                  cell.group.g1=cell.group.g1,
                                  cell.group.g2=cell.group.g2,
                                  psi.method=psi.method,
                                  psi.pval=psi.pval,
                                  psi.delta=psi.delta,
                                  method.de.gene=method.de.gene,
                                  method.adjust.de.gene=method.adjust.de.gene,
                                  downsample=downsample,
                                  seed=seed,
                                  show.progress=show.progress,
                                  mast.method=mast.method,
                                  mast.ebayes=mast.ebayes
                                  )
        
    }
    
}
