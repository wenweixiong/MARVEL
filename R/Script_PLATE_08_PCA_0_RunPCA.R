#' @title Principle component analysis
#'
#' @description Performs principle component analysis on splicing or gene data. This is a wrapper function for \code{RunPCA.PSI} and \code{RunPCA.Exp}.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{ComputePSI} function.
#' @param cell.group.column Character string. The name of the sample metadata column in which the variables will be used to label the cell groups on the PCA.
#' @param cell.group.order Character string. The order of the variables under the sample metadata column specified in \code{cell.group.column} to appear in the PCA cell group legend.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns} and \code{cell.group.order}. If not specified, default \code{ggplot2} colors will be used.
#' @param sample.ids Character strings. Specific cells to plot.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event or gene for the event or gene, respectively, to be included for analysis.
#' @param min.events.pct Numeric value.  Only applicable when \code{level} set to \code{"splicing"}. The minimum percentage of events expressed in a cell, above which, the cell will be retained for analysis. By default, this option is switched off, i.e., \code{NULL}.
#' @param features Character string. Vector of \code{tran_id} or \code{gene_id} for analysis. Should match \code{tran_id} or \code{gene_id} column of \code{MarvelObject$ValidatedSpliceFeature} or \code{MarvelObject$GeneFeature} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param point.alpha Numeric value. Transparency of the data points on reduced dimension space. Take any values between 0 to 1. The smaller the value, the more transparent the data points will be.
#' @param point.stroke Numeric value. The thickness of the outline of the data points. The larger the value, the thicker the outline of the data points.
#' @param level Character string. Indicate \code{"splicing"}, \code{"gene"}, or, \code{"integrated"} for splicing, gene expression analysis, or combined splicing and gene expression analysis, respectively. For \code{"integrated"}, users should run both \code{"splicing"} and \code{"gene"} prior to running \code{"integrated"}.
#' @param method.impute Character string. Only applicable when \code{level} set to \code{"splicing"}. Indicate the method for imputing missing PSI values (low coverage). \code{"random"} method randomly assigns any values between 0-1. \code{"Bayesian"} method uses the posterior PSI computed from the \code{ComputePSI.Posterior} function. Default is \code{"random"}.
#' @param seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}. Ensures imputed values for NA PSIs are reproducible when \code{method.impute} option set to \code{"random"}. Default value is \code{1}.
#' @param pcs Numeric vector. The two principal components (PCs) to plot. Default is the first two PCs. If a vector of 3 is specified, a 3D scatterplot is returned.
#' @param mode Character string. Specify \code{"pca"} for linear dimension reduction analysis or \code{"umap"} for non-linear dimension reduction analysis. Specify \code{"elbow.plot"} to return eigen values. Default is \code{"pca"}.
#' @param seed.umap Numeric value. Only applicable when \code{mode} set to \code{"umap"}. To sure reproducibility of analysis. Default value is \code{42}.
#' @param npc.umap Numeric value. Only applicable when \code{level} set to \code{"splicing"} or \code{"gene"}. Indicate the first number of principal components to use for UMAP . Default value is \code{30}, i.e., the first 30 PCs.
#' @param n.dim Numeric value. Only applicable when \code{level} set to \code{"integrated"}. Indicate the first number of principal components to use for UMAP . Default value is \code{20}, i.e., the first 20 PCs.
#' @param remove.outliers Logical value. If set to \code{TRUE}, outliers will be removed. Outliers defined as data points beyond 1.5 times the interquartile range (IQR) from the 1st and 99th percentile. Default is \code{FALSE}.
#' @param npc.elbow.plot  Numeric value. Only applicable when \code{mode} set to \code{"elbow.plot"}. Incidate the number of PCs to for elbow plot. Default value is \code{50}.
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

RunPCA <- function(MarvelObject,
                   cell.group.column, cell.group.order=NULL, cell.group.colors=NULL,
                   sample.ids=NULL,
                   min.cells=25, min.pct.events=NULL, features,
                   point.size=0.5, point.alpha=0.75, point.stroke=0.1,
                   method.impute="random", seed=1,
                   level,
                   pcs=c(1,2),
                   mode="pca", seed.umap=42, npc.umap=30,
                   n.dim=20,
                   remove.outliers=FALSE,
                   npc.elbow.plot=50
                   ) {

    
    if(level=="splicing") {
        
        RunPCA.PSI(MarvelObject=MarvelObject,
                   cell.group.column=cell.group.column,
                   cell.group.order=cell.group.order,
                   cell.group.colors=cell.group.colors,
                   sample.ids=sample.ids,
                   min.cells=min.cells,
                   min.pct.events=min.pct.events,
                   features=features,
                   point.size=point.size,
                   point.alpha=point.alpha,
                   point.stroke=point.stroke,
                   method.impute=method.impute,
                   seed=seed,
                   pcs=pcs,
                   mode=mode,
                   seed.umap=seed.umap,
                   npc.umap=npc.umap,
                   remove.outliers=remove.outliers,
                   npc.elbow.plot=npc.elbow.plot
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
                   point.stroke=point.stroke,
                   pcs=pcs,
                   mode=mode,
                   seed.umap=seed.umap,
                   npc.umap=npc.umap,
                   remove.outliers=remove.outliers,
                   npc.elbow.plot=npc.elbow.plot
                   )
         
    } else if(level=="integrated") {
        
        RunPCA.PSI.Exp(MarvelObject=MarvelObject,
                   n.dim=n.dim,
                   cell.group.column=cell.group.column,
                   cell.group.order=cell.group.order,
                   cell.group.colors=cell.group.colors,
                   sample.ids=sample.ids,
                   point.size=point.size,
                   point.alpha=point.alpha,
                   point.stroke=point.stroke,
                   seed=seed.umap
                   )
        
    }
    
}
