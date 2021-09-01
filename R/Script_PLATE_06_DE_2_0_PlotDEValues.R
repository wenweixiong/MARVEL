#' @title Plot Differential Analysis Results
#'
#' @description
#' \code{PlotDEValues} plots the differential splicing or gene expression analysis results.
#'
#' @details
#' This function plot the differential splicing or gene expression analysis results and can annotate the user-defined significant data points and gene names.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues} function.
#' @param de.p.val.adj Numeric value. Adjusted p-value below which the splcing events or genes are considered as statistically significant and will consequently be color-annotated on the plot.
#' @param level Character string. Indicate \code{"splicing.distance"} if the percent spliced-in (PSI) values' distribution was previously tested between 2 groups of cells using the \code{CompareValues} function. Statistical tests for distribution include Kolmogorov-Smirnov, Kuiper, and Anderson-Darling test. Indicate \code{"splicing.mean"} or \code{gene} if the PSI or gene expression values' mean was previously tested between 2 groups of cells using the \code{CompareValues} function. Statistical tests for comparing mean are t-test and Wilcoxon rank-sum test.
#' @param de.mean.diff Numeric value. Only applicable when \code{level} set to \code{"splicing.mean"}. The positive (and negative) value specified above (and below) which the splicing events are considered to be statistically significant and will consequently be color-annotated on the plot.
#' @param log2fc Numeric value. Only applicable when \code{level} set to \code{"gene"}. The positive (and negative) value specified above (and below) which the genes are considered to be statistically significant and will consequently be color-annotated on the plot.
#' @param anno Logical value. If set to \code{TRUE}, the specific gene names will be annotated on the plot. Speficified together with \code{anno.x.pos}, \code{anno.y.pos}, \code{anno.x.neg}, and \code{anno.y.neg}.
#' @param anno.x.pos Numeric value. Only applicable if \code{anno} set to TRUE and when \code{level} set to \code{"splicing.mean"} or \code{gene}. The value above on the x-axis which the gene names will be annotated on the plot. Specified together with \code{anno.y.pos}, \code{anno.x.neg}, and \code{anno.y.neg}.
#' @param anno.y.pos Numeric value. Only applicable if \code{anno} set to TRUE and when \code{level} set to \code{"splicing.mean"} or \code{gene}. The value above on the y-axis which the gene names will be annotated on the plot. Specified together with \code{anno.x.pos}, \code{anno.x.neg}, and \code{anno.y.neg}.
#' @param anno.x.neg Numeric value. Only applicable if \code{anno} set to TRUE and when \code{level} set to \code{"splicing.mean"} or \code{gene}. The value below on the x-axis which the gene names will be annotated on the plot. Specified together with \code{anno.x.pos}, \code{anno.y.pos}, and \code{anno.y.neg}.
#' @param anno.y.neg Numeric value. Only applicable if \code{anno} set to TRUE and when \code{level} set to \code{"splicing.mean"} or \code{gene}. The value above on the y-axis which the gene names will be annotated on the plot. Specified together with \code{anno.y.pos}, \code{anno.x.neg}, and \code{anno.x.pos}.
#' @param n.top Numeric value. Only applicable if \code{anno} set to TRUE and when \code{level} set to \code{"splicing.distance"}. Top n differentially expressed splicing events are annotated on the plot.
#' @param label.size Numeric value. Only applicable if \code{anno} set to TRUE. Size of the gene name labels.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$DE$PSI$Plot} or \code{MarvelObject$DE$Exp$Plot} for splicing or gene data, respectively.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import ggplot2
#' @import ggrepel
#' @import scales
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- PlotDEValues(MarvelObject=marvel,
#'                        de.p.val.adj=0.10,
#'                        de.mean.diff=0.01,
#'                        level="splicing.mean",
#'                        anno=FALSE
#'                        )
#'
#' # Check output
#' marvel$DE$PSI$Plot

PlotDEValues <- function(MarvelObject, de.p.val.adj, level,
                         de.mean.diff=NULL, log2fc=NULL,
                         anno=FALSE, anno.x.pos=NULL, anno.y.pos=NULL, anno.x.neg=NULL, anno.y.neg=NULL, n.top=NULL, label.size=2.5
                         ) {

    
    if(level=="splicing.mean") {
        
        PlotDEValues.PSI.Mean(MarvelObject=MarvelObject,
                              de.p.val.adj=de.p.val.adj,
                              de.mean.diff=de.mean.diff,
                              anno=anno,
                              anno.x.pos=anno.x.pos,
                              anno.y.pos=anno.y.pos,
                              anno.x.neg=anno.x.neg,
                              anno.y.neg=anno.y.neg,
                              label.size=label.size
                              )
                              
    } else if(level=="splicing.distance") {
        
        PlotDEValues.PSI.Distance(MarvelObject=MarvelObject,
                                 de.p.val.adj,
                                 anno=anno,
                                 n.top=n.top,
                                 label.size=label.size
                                 )
        
                
    } else if(level=="gene") {
        
        PlotDEValues.Exp(MarvelObject=MarvelObject,
                         de.p.val.adj=de.p.val.adj,
                         log2fc=log2fc,
                         anno=anno,
                         anno.x.pos=anno.x.pos,
                         anno.y.pos=anno.y.pos,
                         anno.x.neg=anno.x.neg,
                         anno.y.neg=anno.y.neg,
                         label.size=label.size
                         )
        
    }

}
