#' @title Differential gene expression analysis for differentially spliced genes
#'
#' @description Performs differential gene expression analysis between 2 groups of cells only on differentially spliced genes.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2.
#' @param downsample Logical value. If set to \code{TRUE}, the number of cells in each cell group will be downsampled to the sample size of the smaller cell group so that both cell groups will have the sample size prior to differential expression analysis. Default is \code{FALSE}.
#' @param seed Numeric value. The seed number for the random number generator to ensure reproducibility during during down-sampling of cells when \code{downsample} set to \code{TRUE}.
#' @param psi.method Vector of character string(s). To include significant events from these method(s) for differential gene expression analysis.
#' @param psi.pval Vector of numeric value(s). The adjusted p-value, below which, the splicing event is considered differentially spliced, and the corresponding genes will be included for differential gene expression analysis.
#' @param psi.delta Numeric value. The absolute difference in mean PSI values between \code{cell.group.g1} and \code{cell.group.g1}, above which, the splicing event is considered differentially spliced, and the corresponding genes will be included for differential gene expression analysis.
#' @param method.de.gene Character string. Same as \code{method} in \code{CompareValues} function.
#' @param method.adjust.de.gene Character string. Same as \code{method} in \code{CompareValues} function.
#' @param show.progress Logical value. If set to \code{TRUE}, progress bar will be displayed so that users can estimate the time needed for differential analysis. Default value is \code{TRUE}.
#' @param mast.method Character string. As per the \code{method} option of the \code{zlm} function from the \code{MAST} package. Default is \code{"bayesglm"}, other options are \code{"glm"} and \code{"glmer"}.
#' @param mast.ebayes Logical value. As per the \code{ebayes} option of the \code{zlm} function from the \code{MAST} package. Default is \code{TRUE}.
#'
#' @return An object of class S3 new slot \code{MarvelObject$DE$Exp$Table}.
#'
#' @importFrom plyr join
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
#' marvel.demo <- CompareValues.Exp.Spliced(MarvelObject=marvel.demo,
#'                                          cell.group.g1=cell.group.g1,
#'                                          cell.group.g2=cell.group.g2,
#'                                          psi.method="ad",
#'                                          psi.pval=0.10,
#'                                          psi.delta=0,
#'                                          method.de.gene="t.test",
#'                                          method.adjust.de.gene="fdr",
#'                                          show.progress=FALSE
#'                                          )
#'
#' # Check output
#' head(marvel.demo$DE$Exp.Spliced$Table)

CompareValues.Exp.Spliced <- function(MarvelObject, cell.group.g1=NULL, cell.group.g2=NULL, psi.method, psi.pval, psi.delta, method.de.gene="wilcox", method.adjust.de.gene="fdr", downsample=FALSE, seed=1, show.progress=TRUE, mast.method="bayesglm", mast.ebayes=TRUE) {

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
    
    # Example arguments
    #MarvelObject <- marvel
    #cell.group.g1 <- cell.group.g1
    #cell.group.g2 <- cell.group.g2
    #psi.method <- c("ad", "dts")
    #psi.pval <- c(0.1, 0.1)
    #psi.delta <- 0
    #method.de.gene <- "mast"
    #method.adjust.de.gene <- "fdr"
    #downsample <- FALSE
    #show.progress <- TRUE
    #mast.method="bayesglm"
    #mast.ebayes=TRUE
    
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
    
    # Perform DE analysis
    object <- CompareValues(MarvelObject=MarvelObject,
                            cell.group.g1=cell.group.g1,
                            cell.group.g2=cell.group.g2,
                            min.cells=3,
                            method=method.de.gene,
                            method.adjust=method.adjust.de.gene,
                            level="gene",
                            custom.gene_ids=df$gene_id,
                            downsample=downsample,
                            seed=seed,
                            show.progress=show.progress,
                            mast.method="bayesglm",
                            mast.ebayes=TRUE
                            )
                            
    ############################################################
    
    # Save to new slot
    MarvelObject$DE$Exp.Spliced$Table <- object$DE$Exp.Custom$Table
    
    # Return final object
    return(MarvelObject)
        
}
