#' @title Differential Splicing and Gene Expression Analysis
#'
#' @description
#' \code{CompareValues} performs differentially splicing and gene expression analysis between 2 groups of cells.
#'
#' @details
#' This function compares the percent spliced-in (PSI) and gene expression values between 2 groups of cells.
#'
#' @param MarvelObject S3 object generated from \code{ComputePSI} function.
#' @param cell.type.columns.1 Character string. To indicate which columns in the \code{MarvelObject$SplicePheno} or  \code{MarvelObject$GenePheno} slot to refer to when filtering samples (cells) for analysis for splicing or gene data, respectively. This is for the first group of cells (reference group).
#' @param cell.type.variables.1 List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{cell.type.columns} argument. This is for the first group of cells (reference group).
#' @param cell.type.columns.2 Character string. To indicate which columns in the \code{MarvelObject$SplicePheno} or  \code{MarvelObject$GenePheno} slot to refer to when filtering samples (cells) for analysis for splicing or gene data, respectively. This is for the 2nd group of cells (non-reference group).
#' @param cell.type.variables.2 List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{cell.type.columns} argument. This is for the 2nd group of cells  (non-reference group).
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event for the event to be included for differential splicing analysis.
#' @param method Character string. Statistical test to compare the 2 groups of cells. \code{"ks"}, \code{"kuiper"}, \code{"ad"}, \code{"wilcox"}, and \code{"t.test"} for Kolmogorov-Smirnov, Kuiper, Anderson-Darling, Wilcox, and t-test, respectively.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for differential splicing or gene expression analysis, respectively
#' @param event.type Character string. Only applicable when \code{level} set to \code{"splicing"}. Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot. When \code{level} set to \code{"splicing"} or \code{"gene"}, results are returned to \code{$DE$PSI} or \code{$DE$Exp} slot, respectively.
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
#' marvel <- CompareValues(MarvelObject=marvel,
#'                         cell.type.columns.1=c("cell.type"),
#'                         cell.type.variables.1=list("iPSC"),
#'                         cell.type.columns.2=c("cell.type"),
#'                         cell.type.variables.2=list("Endoderm"),
#'                         n.cells=2,
#'                         method="t.test",
#'                         method.adjust="fdr",
#'                         level="splicing",
#'                         event.type=c("SE", "MXE", "RI", "A5SS", "A3SS")
#'                         )
#'
#' # Check output
#' marvel$DE$PSI$Table[1:5, ]

CompareValues <- function(MarvelObject, cell.type.columns.1, cell.type.variables.1, cell.type.columns.2, cell.type.variables.2, n.cells, method, method.adjust, level, event.type=NULL) {

    # Define arguments
    MarvelObject <- MarvelObject
    cell.type.columns.1 <- cell.type.columns.1
    cell.type.variables.1 <- cell.type.variables.1
    cell.type.columns.2 <- cell.type.columns.2
    cell.type.variables.2 <- cell.type.variables.2
    n.cells <- n.cells
    method <- method
    method.adjust <- method.adjust
    level <- level
    event.type <- event.type
    
    if(level=="splicing") {
        
        CompareValues.PSI(MarvelObject=MarvelObject,
                          cell.type.columns.1=cell.type.columns.1,
                          cell.type.variables.1=cell.type.variables.1,
                          cell.type.columns.2=cell.type.columns.2,
                          cell.type.variables.2=cell.type.variables.2,
                          n.cells=n.cells,
                          method=method,
                          method.adjust=method.adjust,
                          event.type=event.type
                          )

    } else if(level=="gene") {
        
        CompareValues.Exp(MarvelObject=MarvelObject,
                          cell.type.columns.1=cell.type.columns.1,
                          cell.type.variables.1=cell.type.variables.1,
                          cell.type.columns.2=cell.type.columns.2,
                          cell.type.variables.2=cell.type.variables.2,
                          n.cells=n.cells,
                          method=method,
                          method.adjust=method.adjust
                          )
        
    }
    
}
