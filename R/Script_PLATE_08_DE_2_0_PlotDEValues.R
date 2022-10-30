#' @title Plot differential splicing and gene expression analysis results
#'
#' @description Volcano plot of differential splicing and gene expression analysis results. This is a wrapper function for \code{PlotDEValues.PSI.Mean}, \code{PlotDEValues.Exp.Global}, and \code{PlotDEValues.Exp.Spliced}.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues} function.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param pval Numeric value. Only applicable when \code{level} set to \code{"splicing.mean"}, \code{"splicing.distance"}, and \code{"gene.global"}. Adjusted p-value below which the splcing events or genes are considered as statistically significant and will consequently be color-annotated on the plot.
#' @param level Character string. Indicate \code{"splicing.distance"} if the percent spliced-in (PSI) values' distribution was previously tested between 2 groups of cells using the \code{CompareValues} function. Statistical tests for distribution include Kolmogorov-Smirnov, Kuiper, and Anderson-Darling test. Indicate \code{"splicing.mean"} or \code{gene} if the PSI or gene expression values' mean was previously tested between 2 groups of cells using the \code{CompareValues} function. Statistical tests for comparing mean are t-test and Wilcoxon rank-sum test.
#' @param delta Numeric value. Only applicable when \code{level} set to \code{"splicing.mean"}. The positive (and negative) value specified above (and below) which the splicing events are considered to be statistically significant and will consequently be color-annotated on the plot.
#' @param log2fc Numeric value. Only applicable when \code{level} set to \code{"gene.global"}. The positive (and negative) value specified above (and below) which the genes are considered to be statistically significant and will consequently be color-annotated on the plot.
#' @param psi.pval Numeric value. Only applicable when \code{level} set to \code{"gene.spliced"}. The adjusted p-value from differential splicing analysis, below which, the splicing event is considered differentially spliced. Default is \code{0.1}.
#' @param psi.delta Numeric value. Only applicable when \code{level} set to \code{"gene.spliced"}. The absolute differences in average PSI value between two cell groups from differential splicing analysis, above which, the splicing event is considered differentially spliced.  Default is \code{0}.
#' @param gene.pval Numeric value. Only applicable when \code{level} set to \code{"gene.spliced"}. The adjusted p-value from differential gene expression analysis, below which, the gene is considered differentially expressed. Default is \code{0.1}.
#' @param gene.log2fc Numeric value. Only applicable when \code{level} set to \code{"gene.spliced"}. The absolute log2 fold change in gene expression betwene two cell groups from differential splicing analysis, above which, the gene is considered differentially expressed. Default is \code{0.5}.
#' @param point.size Numeric value. Size of data points. Default is \code{1}.
#' @param anno Logical value. If set to \code{TRUE}, the specific gene names or splicing events will be annotated on the plot.
#' @param anno.tran_id Vector of character strings. When \code{anno} set to \code{TRUE}, the coordinates of the splicing events to be annotated on the plot.
#' @param anno.gene_short_name Vector of character strings. When \code{anno} set to \code{TRUE}, the gene names to be annotated on the plot.
#' @param label.size Numeric value. Only applicable if \code{anno} set to TRUE. Size of the gene name labels.
#' @param y.upper.offset Numeric value. The value in -log10(p-value) to increase the upper limit of the y-axis. To be used when \code{anno} set to TRUE so that gene labels will not be truncated at the upper limit of the y-axis.
#' @param xlabel.size Numeric value. Font size of the xtick labels. Default is \code{8}.
#' @param point.alpha Numeric value. Only applicable when \code{level} set to \code{"splicing.mean.g2vsg1"}. Transpancy of data points. Default is \code{1}.
#' @param event.types Vector of character string(s). Only applicable when \code{level} set to \code{"splicing.mean.g2vsg1"}. The specific splicing event to plot. May take any one or more of the following values \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, \code{"A3SS"}, \code{"AFE"}, and  \code{"ALE"}.
#' @param event.types.colors Vector of character string(s). Only applicable when \code{level} set to \code{"splicing.mean.g2vsg1"}. Customise colors as per splicing event type specified in \code{event.types} option. Should be of same length as \code{event.types} option.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$DE$PSI$Plot[["method"]]} when level set to \code{"splicing.mean"} or \code{"splicing.distance"} or \code{MarvelObject$DE$Exp.Global$Table} and \code{MarvelObject$DE$Exp.Global$Plot} when level set to \code{"gene.global"} or \code{MarvelObject$DE$Exp.Spliced$Table} and \code{MarvelObject$DE$Exp.Spliced$Plot} when level set to \code{"gene.spliced"}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @import scales
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- PlotDEValues(MarvelObject=marvel.demo,
#'                             method="ad",
#'                             pval=0.10,
#'                             level="splicing.distance"
#'                             )
#'
#' # Check output
#' marvel.demo$DE$PSI$Plot[["ad"]]

PlotDEValues <- function(MarvelObject, method=NULL, pval, level,
                         delta=NULL, log2fc=NULL,
                         psi.pval=NULL, psi.delta=NULL, gene.pval=NULL, gene.log2fc=NULL,
                         point.size=1, xlabel.size=8, point.alpha=1,
                         anno=FALSE, anno.gene_short_name=NULL, anno.tran_id=NULL, label.size=2.5,
                         y.upper.offset=5,
                         event.types=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"), event.types.colors=NULL
                         ) {

    
    if(level=="splicing.mean") {
        
        PlotDEValues.PSI.Mean(MarvelObject=MarvelObject,
                              method=method,
                              pval=pval,
                              delta=delta,
                              point.size=point.size,
                              xlabel.size=xlabel.size,
                              anno=anno,
                              anno.tran_id=anno.tran_id,
                              label.size=label.size,
                              y.upper.offset=y.upper.offset
                              )
                              
    } else if(level=="splicing.distance") {
        
        PlotDEValues.PSI.Distance(MarvelObject=MarvelObject,
                                  method=method,
                                  pval,
                                  point.size=point.size,
                                  xlabel.size=xlabel.size,
                                  anno=anno,
                                  anno.tran_id=anno.tran_id,
                                  label.size=label.size,
                                  y.upper.offset=y.upper.offset
                                  )
        
                
    } else if(level=="gene.global") {
        
        PlotDEValues.Exp.Global(MarvelObject=MarvelObject,
                                pval=pval,
                                log2fc=log2fc,
                                point.size=point.size,
                                xlabel.size=xlabel.size,
                                anno=anno,
                                anno.gene_short_name=anno.gene_short_name,
                                label.size=label.size,
                                y.upper.offset=y.upper.offset
                                )
        
    } else if(level=="gene.spliced") {
        
        PlotDEValues.Exp.Spliced(MarvelObject=MarvelObject,
                                 method=method,
                                 psi.pval=psi.pval,
                                 psi.delta=psi.delta,
                                 gene.pval=gene.pval,
                                 gene.log2fc=gene.log2fc,
                                 point.size=point.size,
                                 xlabel.size=xlabel.size,
                                 anno=anno,
                                 anno.gene_short_name=anno.gene_short_name,
                                 label.size=label.size,
                                 y.upper.offset=y.upper.offset
                                 )
                                
    } else if(level=="splicing.mean.g2vsg1") {
        
        PlotDEValues.PSI.Mean.g2vsg1(MarvelObject=MarvelObject,
                                     method=method,
                                     pval=pval,
                                     delta=delta,
                                     point.size=point.size,
                                     xlabel.size=xlabel.size,
                                     anno=anno,
                                     anno.tran_id=anno.tran_id,
                                     label.size=label.size,
                                     event.types=event.types,
                                     point.alpha=point.alpha,
                                     event.types.colors=event.types.colors
                                     )
        
        
        
    }

}
