#' @title Principle component analysis
#'
#' @description Performs principle component analysis on splicing or gene data. This is a wrapper function for \code{RunPCA.PSI} and \code{RunPCA.Exp}.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.column Character string. The name of the sample metadata column in which the variables will be used to label the cell groups on the PCA.
#' @param cell.group.order Character string. The order of the variables under the sample metadata column specified in \code{cell.group.column} to appear in the PCA cell group legend.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns} and \code{cell.group.order}. If not specified, default \code{ggplot2} colors will be used.
#' @param sample.ids Character strings. Specific cells to plot.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event or gene for the event or gene, respectively, to be included for analysis.
#' @param features Character string. Vector of \code{tran_id} or \code{gene_id} for analysis. Should match \code{tran_id} or \code{gene_id} column of \code{MarvelObject$ValidatedSpliceFeature} or \code{MarvelObject$GeneFeature} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param point.alpha Numeric value. Transparency of the data points on reduced dimension space. Take any values between 0 to 1. The smaller the value, the more transparent the data points will be.
#' @param point.stroke Numeric value. The thickness of the outline of the data points. The larger the value, the thicker the outline of the data points.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for splicing or gene expression analysis, respectively
#' @param seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}.  Ensures imputed values for NA PSIs are reproducible.
#' @param method.impute Character string. Only applicable when \code{level} set to \code{"splicing"}. Indicate the method for imputing missing PSI values (low coverage). \code{"random"} method randomly assigns any values between 0-1. \code{"population.mean"} method uses the mean PSI value for each cell population. Default option is \code{"population.mean"}.
#' @param cell.group.column.impute Character string. Only applicable when \code{method.impute} set to \code{"population.mean"}. The name of the sample metadata column in which the variables will be used to impute missing values.
#'
#' @export
#'
#' @return An object of class S3 with new slots \code{MarvelObject$PCA$PSI$Results}, \code{MarvelObject$PCA$PSI$Plot}, and \code{MarvelObject$PCA$PSI$Plot.Elbow} or \code{MarvelObject$PCA$Exp$Results}, \code{MarvelObject$PCA$Exp$Plot}, and \code{MarvelObject$PCA$Exp$Plot.Elbow}, when \code{level} option specified as \code{"splicing"} or \code{"gene"}, respectively.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' # Define splicing events for analysis
#' df <- do.call(rbind.data.frame, marvel.demo$PSI)
#' tran_ids <- df$tran_id
#'
#' # PCA
#' marvel.demo <- RunPCA(MarvelObject=marvel.demo,
#'                       sample.ids=marvel.demo$SplicePheno$sample.id,
#'                       cell.group.column="cell.type",
#'                       cell.group.order=c("iPSC", "Endoderm"),
#'                       cell.group.colors=NULL,
#'                       min.cells=5,
#'                       features=tran_ids,
#'                       level="splicing",
#'                       point.size=2
#'                       )
#' 
#' # Check outputs
#' head(marvel.demo$PCA$PSI$Results$ind$coord)
#' marvel.demo$PCA$PSI$Plot

RunPCA <- function(MarvelObject, cell.group.column, cell.group.order=NULL, cell.group.colors=NULL,
                   sample.ids=NULL,
                   min.cells=25, features,
                   point.size=0.5, point.alpha=0.75, point.stroke=0.1,
                   seed=1, method.impute="random", cell.group.column.impute=NULL,
                   level
                   ) {

    
    if(level=="splicing") {
        
        RunPCA.PSI(MarvelObject=MarvelObject,
                   cell.group.column=cell.group.column,
                   cell.group.order=cell.group.order,
                   cell.group.colors=cell.group.colors,
                   sample.ids=sample.ids,
                   min.cells=min.cells,
                   features=features,
                   point.size=point.size,
                   point.alpha=point.alpha,
                   point.stroke=point.stroke,
                   seed=seed,
                   method.impute=method.impute,
                   cell.group.column.impute=cell.group.column.impute
                   )

    } else if(level=="gene") {
        
        RunPCA.Exp(MarvelObject=MarvelObject,
                   cell.group.column=cell.group.column,
                   cell.group.order=cell.group.order,
                   cell.group.colors=cell.group.colors,
                   sample.ids=sample.ids,
                   min.cells=min.cells,
                   features=features,
                   point.size=point.size,
                   point.alpha=point.alpha,
                   point.stroke=point.stroke
                   )
        
    }
    
}
