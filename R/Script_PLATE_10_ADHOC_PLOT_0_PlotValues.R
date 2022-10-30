#' @title Plot percent spliced-in (PSI) or gene expression values
#'
#' @description Plots percent spliced-in (PSI) or gene expression values across different groups of cells. This is a wrapper function for \code{PlotValues.Exp} and \code{PlotValues.PSI}.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.list List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group. The name of the element will be the cell group label.
#' @param feature Character string. \code{tran_id} or \code{gene_id} for plotting. Should match \code{tran_id} or \code{gene_id} column of \code{MarvelObject$ValidatedSpliceFeature} or \code{MarvelObject$GeneFeature} slot when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param maintitle Character string. Column to use as plot main title as per \code{MarvelObject$ValidatedSpliceFeature} or \code{MarvelObject$GeneFeature} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively. Default is \code{"gene_short_name"} column.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function. Default is 8.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for PSI or gene expression value plotting, respectively.
#' @param min.cells Numeric value. Only applicable when \code{level} set to \code{"splicing"}. The minimum no. of cells expressing the splicing event to be included for analysis.
#' @param sigma.sq Numeric value. Only applicable when \code{level} set to \code{"splicing"}. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality. Please refer to \code{AssignModality} function help page for more details. Default is 0.001.
#' @param bimodal.adjust Logical. Only applicable when \code{level} set to \code{"splicing"}. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality. Please refer to \code{AssignModality} function help page for more details.
#' @param seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.column Character string. Only applicable when \code{level} set to \code{"splicing"}. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details. Default is \code{"modality.bimodal.adj"}.
#' @param scale.y.log Logical value. Only applicable when \code{level} set to \code{"splicing"}. If set to \code{TRUE}, the y-axis of will log10-scaled. Useful when most PSI values are extremely small (< 0.02) or big (> 0.98). Default is \code{FALSE}.
#' @param max.cells.jitter Numeric value. Only applicable when \code{level} set to \code{"splicing"}. Maximum number of cells for jitter points. Cells are randomly downsampled to show on jitter plot. Useful when there are large number of cells so that individual jitter points do not overcrowd the violin plot. Specified together with \code{max.cells.jitter.seed}. To disable this option, specify a value large than the number of cells in each cell group.
#' @param max.cells.jitter.seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}. Cells downsampled are reproducible. Specified together with \code{max.cells.jitter}.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns}, \code{cell.type.variable}, and \code{cell.type.labels}. If not specified, default \code{ggplot2} colors will be used.
#' @param point.alpha Numeric value. Transparency of the data points. Takes any values between 0-1. Default value is \code{0.2}.
#'
#' @return An object of class S3 with new slot \code{$adhocPlot$PSI} or \code{MarvelObject$adhocPlot$Exp} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#'
#' @importFrom plyr join
#' @import methods
#' @import ggplot2
#' @import scales
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' # Define cell groups to plot
#' df.pheno <- marvel.demo$SplicePheno
#' cell.group.g1 <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
#' cell.group.g2 <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]
#' cell.group.list <- list(cell.group.g1, cell.group.g2)
#' names(cell.group.list) <- c("iPSC", "Endoderm")
#'
#' # Plot
#' marvel.demo <- PlotValues(MarvelObject=marvel.demo,
#'                           cell.group.list=cell.group.list,
#'                           feature="chr17:8383254:8382781|8383157:-@chr17:8382143:8382315",
#'                           level="splicing",
#'                           min.cells=5,
#'                           xlabels.size=5
#'                           )
#'
#' # Check output
#' marvel.demo$adhocPlot$PSI

PlotValues <- function(MarvelObject, cell.group.list, feature, maintitle="gene_short_name", xlabels.size=8, level, min.cells=NULL, sigma.sq=0.001, bimodal.adjust=NULL, seed=NULL, modality.column="modality.bimodal.adj", scale.y.log=FALSE, max.cells.jitter=10000, max.cells.jitter.seed=1, cell.group.colors=NULL, point.alpha=0.2) {
    
    if(level=="gene") {
        
        PlotValues.Exp(MarvelObject=MarvelObject,
                       cell.group.list=cell.group.list,
                       feature=feature,
                       maintitle=maintitle,
                       xlabels.size=xlabels.size,
                       cell.group.colors=cell.group.colors,
                       point.alpha=point.alpha
                   )

    } else if(level=="splicing") {
        
        PlotValues.PSI(MarvelObject=MarvelObject,
                       cell.group.list=cell.group.list,
                       feature=feature,
                       maintitle=maintitle,
                       xlabels.size=xlabels.size,
                       max.cells.jitter=max.cells.jitter,
                       max.cells.jitter.seed=max.cells.jitter.seed,
                       min.cells=min.cells,
                       sigma.sq=sigma.sq,
                       bimodal.adjust=bimodal.adjust,
                       seed=seed,
                       modality.column=modality.column,
                       scale.y.log=scale.y.log,
                       cell.group.colors=cell.group.colors,
                       point.alpha=point.alpha
                   )
        
    }

    
}
