#' @title Plot Percent Spliced-in (PSI) and Gene Expression Values
#'
#' @description
#' \code{PlotValues} plots percent spliced-in (PSI) and gene expression values across different groups of cells.
#'
#' @details
#' This function plots percent spliced-in (PSI) and gene expression values across different groups of cells. Violin plot is used for PSI values while boxplot is used for gene expression values.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.types Character string. To indicate which groups of cells that will be used for plotting. Group names should match those in \code{cell.type} column of \code{$SplicePheno} or \code{GenePheno} slot for splicing or gene expression data, respectively.
#' @param feature Character string. \code{tran_id} or \code{gene_id} for plotting. Should match \code{tran_id} or \code{gene_id} column of \code{$ValidatedSpliceFeature} or \code{GeneFeature} slot when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param maintitle Character string. Column to use as plot main title as per \code{ValidatedSpliceFeature} or \code{GeneFeature} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for PSI or gene expression value plotting, respectively.
#' @param n.cells Numeric value. Only applicable when \code{level} set to \code{"splicing"}. The minimum no. of cells expressing the splicing event to be included for analysis.
#' @param sigma.sq Numeric value. Only applicable when \code{level} set to \code{"splicing"}. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality. Please refer to \code{AssignModality} function help page for more details.
#' @param bimodal.adjust Logical. Only applicable when \code{level} set to \code{"splicing"}. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality. Please refer to \code{AssignModality} function help page for more details.
#' @param seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.column Character string. Only applicable when \code{level} set to \code{"splicing"}. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{$adhocPlot$PSI} or \code{$adhocPlot$Gene} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @examples
#'
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' feature <- marvel$SpliceFeature$SE$tran_id[1]
#'
#' marvel <- PlotValues(MarvelObject=marvel,
#'                      cell.types=c("iPSC", "Endoderm"),
#'                      feature=feature,
#'                      maintitle="gene_short_name",
#'                      xlabels.size=12,
#'                      level="splicing",
#'                      n.cells=25,
#'                      sigma.sq=0.001,
#'                      bimodal.adjust=TRUE,
#'                      seed=1,
#'                      modality.column="modality.bimodal.adj"
#'                      )
#'
#' marvel$adhocPlot$PSI

PlotValues <- function(MarvelObject, cell.types, feature, maintitle, xlabels.size=NULL, level, n.cells=NULL, sigma.sq=NULL, bimodal.adjust=NULL, seed=NULL, modality.column=NULL) {
    
    if(level=="gene") {
        
        PlotValues.Exp(MarvelObject=MarvelObject,
                       cell.types=cell.types,
                       feature=feature,
                       maintitle=maintitle
                   )

    } else if(level=="splicing") {
        
        PlotValues.PSI(MarvelObject=MarvelObject,
                       cell.types=cell.types,
                       feature=feature,
                       maintitle=maintitle,
                       xlabels.size=xlabels.size,
                       n.cells=n.cells,
                       sigma.sq=sigma.sq,
                       bimodal.adjust=bimodal.adjust,
                       seed=seed,
                       modality.column=modality.column
                   )
        
    }

    
}
