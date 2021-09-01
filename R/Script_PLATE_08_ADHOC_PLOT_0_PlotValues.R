#' @title Plot Percent Spliced-in (PSI) or Gene Expression Values
#'
#' @description
#' \code{PlotValues} plots percent spliced-in (PSI) or gene expression values across different groups of cells.
#'
#' @details
#' This function plots percent spliced-in (PSI) or gene expression values across different groups of cells. Violin plot is used for PSI values while boxplot is used for gene expression values.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.type.columns List. To indicate which columns in the \code{MarvelObject$SplicePheno} slot to refer to when filtering samples (cells) for analysis. Each item in the list defines the columns for each group to plot.
#' @param cell.type.variables List. To indicate which specific variables of the corresponding columns to keep the samples (cells). Each item in the list should be the same length as the corresponding item in the \code{cell.type.columns} list.
#' @param cell.type.labels Character string. To indicate the cell group labels on the x-axis. Should be same length as the number of items in \code{cell.type.columns} and \code{cell.type.variables} lists.
#' @param feature Character string. \code{tran_id} or \code{gene_id} for plotting. Should match \code{tran_id} or \code{gene_id} column of \code{MarvelObject$ValidatedSpliceFeature} or \code{MarvelObject$GeneFeature} slot when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param maintitle Character string. Column to use as plot main title as per \code{MarvelObject$ValidatedSpliceFeature} or \code{MarvelObject$GeneFeature} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively. Default is \code{"gene_short_name"} column.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function. Default is 8.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for PSI or gene expression value plotting, respectively.
#' @param n.cells Numeric value. Only applicable when \code{level} set to \code{"splicing"}. The minimum no. of cells expressing the splicing event to be included for analysis.
#' @param sigma.sq Numeric value. Only applicable when \code{level} set to \code{"splicing"}. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality. Please refer to \code{AssignModality} function help page for more details. Default is 0.001.
#' @param bimodal.adjust Logical. Only applicable when \code{level} set to \code{"splicing"}. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality. Please refer to \code{AssignModality} function help page for more details.
#' @param seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.column Character string. Only applicable when \code{level} set to \code{"splicing"}. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details. Default is \code{"modality.bimodal.adj"}.
#' @param scale.y Logical value. Only applicable when \code{level} set to \code{"splicing"}. If set to \code{FALSE}, the y-axis of will not be scaled to 0 to 1, but will instead follow the range of the data values. Useful when most PSI values are extremely small (< 0.02) or big (> 0.98). Default is \code{TRUE}.
#' @param n.cells.jitter.threshold Numeric value. Only applicable when \code{level} set to \code{"splicing"}. Maximum number of cells for jitter points. Cells are randomly downsampled to show on jitter plot. Useful when there are large number of cells so that individual jitter points do not overcrowd the violin plot. Specified together with \code{n.cells.jitter.seed}
#' @param n.cells.jitter.seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}. Cells downsampled are reproducible. Specified together with \code{n.cells.jitter.threshold}.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{$adhocPlot$PSI} or \code{MarvelObject$adhocPlot$Exp} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Define event to plot
#' tran_id <- marvel$DE$PSI$Table$tran_id[1]
#'
#' # Run example
#' marvel <- PlotValues(MarvelObject=marvel,
#'                      cell.type.columns=list(c("cell.type"), c("cell.type")),
#'                      cell.type.variables=list(list("iPSC"), list("Endoderm")),
#'                      cell.type.labels=c("iPSC", "Endoderm"),
#'                      feature=tran_id,
#'                      xlabels.size=8,
#'                      level="splicing",
#'                      n.cells=2,
#'                      bimodal.adjust=TRUE,
#'                      seed=1,
#'                      )
#'
#' # Check output
#' marvel$adhocPlot$PSI

PlotValues <- function(MarvelObject, cell.type.columns, cell.type.variables, cell.type.labels, feature, maintitle="gene_short_name", xlabels.size=8, level, n.cells=NULL, sigma.sq=0.001, bimodal.adjust=NULL, seed=NULL, modality.column="modality.bimodal.adj", scale.y=TRUE, n.cells.jitter.threshold=NULL, n.cells.jitter.seed=NULL) {
    
    if(level=="gene") {
        
        PlotValues.Exp(MarvelObject=MarvelObject,
                       cell.type.columns=cell.type.columns,
                       cell.type.variables=cell.type.variables,
                       cell.type.labels=cell.type.labels,
                       feature=feature,
                       maintitle=maintitle,
                       xlabels.size=xlabels.size
                   )

    } else if(level=="splicing") {
        
        PlotValues.PSI(MarvelObject=MarvelObject,
                       cell.type.columns=cell.type.columns,
                       cell.type.variables=cell.type.variables,
                       cell.type.labels=cell.type.labels,
                       feature=feature,
                       maintitle=maintitle,
                       xlabels.size=xlabels.size,
                       n.cells.jitter.threshold=n.cells.jitter.threshold,
                       n.cells.jitter.seed=n.cells.jitter.seed,
                       n.cells=n.cells,
                       sigma.sq=sigma.sq,
                       bimodal.adjust=bimodal.adjust,
                       seed=seed,
                       modality.column=modality.column,
                       scale.y=scale.y
                   )
        
    }

    
}
