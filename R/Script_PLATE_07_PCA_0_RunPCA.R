#' @title Principle component analysis
#'
#' @description Performs principle component analysis on splicing or gene data. This is a wrapper function for \code{RunPCA.PSI} and \code{RunPCA.Exp}.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.list List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group. The name of the element will be the cell group label.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event or gene for the event or gene, respectively, to be included for analysis.
#' @param features Character string. Vector of \code{tran_id} or \code{gene_id} for analysis. Should match \code{tran_id} or \code{gene_id} column of \code{MarvelObject$ValidatedSpliceFeature} or \code{MarvelObject$GeneFeature} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param point.alpha Numeric value. Transparency of the data points on reduced dimension space. Take any values between 0 to 1. The smaller the value, the more transparent the data points will be.
#' @param point.stroke Numeric value. The thickness of the outline of the data points. The larger the value, the thicker the outline of the data points.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for splicing or gene expression analysis, respectively
#' @param event.type Character string. Only applicable when \code{level} set to \code{"splicing"}.  Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}.  Ensures imputed values for NA PSIs are reproducible.
#' @param method.impute Character string. Only applicable when \code{level} set to \code{"splicing"}. Indicate the method for imputing missing PSI values (low coverage). \code{"random"} method randomly assigns any values between 0-1. \code{"population.mean"} method uses the mean PSI value for each cell population. Default option is \code{"population.mean"}.
#' @param retrieve.non.outliers Logical. If set to \code{TRUE}, this function will retrieve \code{sample.id} of non-outliers based on the intial PCA. Define the non-outliers based on the initial PCA coordinates. Use in conjunction with arguments \code{pc1.min}, \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc2.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc1.max}, and \code{pc2.max}.
#' @param pc2.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.min}.
#' @param remove.outliers Logical. If set to \code{TRUE}, re-run PCA by only including non-outliers. Use after running the function with \code{retrieve.non.outliers} set to \code{TRUE}.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns}, \code{cell.type.variable}, and \code{cell.type.labels}. If not specified, default \code{ggplot2} colors will be used.
#'
#' @export
#'
#' @return An object of class S3 with new slots \code{MarvelObject$PCA$PSI$Results}, \code{MarvelObject$PCA$PSI$Plot}, and \code{MarvelObject$PCA$PSI$Plot.Elbow} or \code{MarvelObject$PCA$Exp$Results}, \code{MarvelObject$PCA$Exp$Plot}, and \code{MarvelObject$PCA$Exp$Plot.Elbow}, when \code{level} option specified as \code{"splicing"} or \code{"gene"}, respectively.
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import FactoMineR
#' @import factoextra
#'
#' @export

RunPCA <- function(MarvelObject, cell.group.list, min.cells=25, features,
                   point.size=0.5, point.alpha=0.75, point.stroke=0.1, level, event.type=NULL,
                   seed=1, method.impute="population.mean",
                   retrieve.non.outliers=FALSE, pc1.min=NULL, pc1.max=NULL, pc2.min=NULL, pc2.max=NULL,
                   remove.outliers=FALSE, cell.group.colors=NULL
                   ) {

    
    if(level=="splicing") {
        
        RunPCA.PSI(MarvelObject=MarvelObject,
                   cell.group.list=cell.group.list,
                   min.cells=min.cells,
                   features=features,
                   point.size=point.size,
                   point.alpha=point.alpha,
                   point.stroke=point.stroke,
                   event.type=event.type,
                   seed=seed,
                   method.impute=method.impute,
                   retrieve.non.outliers=retrieve.non.outliers,
                   pc1.min=pc1.min,
                   pc1.max=pc1.max,
                   pc2.min=pc2.min,
                   pc2.max=pc2.max,
                   remove.outliers=remove.outliers,
                   cell.group.colors=cell.group.colors
                   )

    } else if(level=="gene") {
        
        RunPCA.Exp(MarvelObject=MarvelObject,
                   cell.group.list=cell.group.list,
                   min.cells=min.cells,
                   features=features,
                   point.size=point.size,
                   point.alpha=point.alpha,
                   point.stroke=point.stroke,
                   retrieve.non.outliers=retrieve.non.outliers,
                   pc1.min=pc1.min,
                   pc1.max=pc1.max,
                   pc2.min=pc2.min,
                   pc2.max=pc2.max,
                   remove.outliers=remove.outliers,
                   cell.group.colors=cell.group.colors
                   )
        
    }
    
}
