#' @title Differential Splicing and Gene Expression Analysis
#'
#' @description
#' \code{CompareValues} performs differentially splicing and gene expression analysis between 2 groups of cells.
#'
#' @details
#' This function compares the percent spliced-in (PSI) and gene expression values between 2 groups of cells.
#'
#' @param MarvelObject S3 object generated from \code{ComputePSI} function.
#' @param cell.types Character string. To indicate which 2 groups of cells that will be used for differential splicing analysis. Group names should match those in \code{cell.type} column of \code{$SplicePheno} slot.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event for the event to be included for differential splicing analysis.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Wilcox, and t-test, respectively. We advice \code{"ks"} for PSI comparison while \code{"wilcox"} or \code{"t.test"} for gene expression comparison.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for differential splicing or gene expression analysis, respectively
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot. When \code{level} set to \code{"splicing"} or \code{"gene"}, results are returned to \code{$DE$PSI} or \code{$DE$Exp} slot, respectively.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' marvel <- CompareValues(MarvelObject=marvel,
#'                         cell.types=c("iPSC", "Endoderm"),
#'                         n.cells=25,
#'                         method="t.test",
#'                         method.adjust="fdr",
#'                         level="splicing")
#'
#' marvel$DE$PSI[1:5, ]

CompareValues <- function(MarvelObject, cell.types, n.cells, method, method.adjust, level) {

    # Define arguments
    MarvelObject <- MarvelObject
    cell.types <- cell.types
    n.cells <- n.cells
    method <- method
    method.adjust <- method.adjust
    level <- level
    
    if(level=="splicing") {
        
        CompareValues.PSI(MarvelObject=MarvelObject,
                   cell.types=cell.types,
                   n.cells=n.cells,
                   method=method,
                   method.adjust=method.adjust
                   )

    } else if(level=="gene") {
        
        CompareValues.Exp(MarvelObject=MarvelObject,
                   cell.types=cell.types,
                   n.cells=n.cells,
                   method=method,
                   method.adjust=method.adjust
                   )
        
    }
    
}
