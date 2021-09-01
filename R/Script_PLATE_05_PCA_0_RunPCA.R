#' @title Principle Component Analysis
#'
#' @description
#' \code{RunPCA} performs principle component analysis on splicing or gene data.
#'
#' @details
#' This function performs principle component analysis on splicing or gene data and visualise cells on a reducted dimension space, i.e. 2D.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.type.columns Character string. To indicate which columns in the \code{MarvelObject$SplicePheno} or  \code{MarvelObject$GenePheno} slot to refer to when filtering samples (cells) for analysis for splicing or gene data, respectively.
#' @param cell.type.variables List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{cell.type.columns} argument.
#' @param cell.type.labels Character string. To indicate the cell group labels. Should be same length as the number of items in \code{cell.type.columns} and \code{cell.type.variables} lists.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event or gene for the event or gene, respectively, to be included for analysis.
#' @param features Character string. Vector of \code{tran_id} or \code{gene_id} for analysis. Should match \code{tran_id} or \code{gene_id} column of \code{MarvelObject$ValidatedSpliceFeature} or \code{MarvelObject$GeneFeature} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for splicing or gene expression analysis, respectively
#' @param event.type Character string. Only applicable when \code{level} set to \code{"splicing"}.  Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}.  Ensures imputed values for NA PSIs are reproducible.
#' @param retrieve.non.outliers Logical. If set to \code{TRUE}, this function will retrieve \code{sample.id} of non-outliers based on the intial PCA. Define the non-outliers based on the initial PCA coordinates. Use in conjunction with arguments \code{pc1.min}, \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc2.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc1.max}, and \code{pc2.max}.
#' @param pc2.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.min}.
#' @param remove.outliers Logical. If set to \code{TRUE}, re-run PCA by only including non-outliers. Use after running the function with \code{retrieve.non.outliers} set to \code{TRUE}.

#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$PCA$PSI} or \code{MarvelObject$PCA$Gene} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively. Contains 2D scatterplot in \code{MarvelObject$PCA$PSI$Plot} or \code{MarvelObject$PCA$Gene$Plot} and the corresponding x- and y-coordinates for each sample in \code{MarvelObject$PCA$PSI$Results} or \code{MarvelObject$PCA$Gene$Results}, respectively.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import FactoMineR
#' @import factoextra
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Define events to reduce dimension
#' . <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
#' tran_ids <- .$tran_id
#'
#' # Run example
#' marvel <- RunPCA(MarvelObject=marvel,
#'                  cell.type.columns=list(c("cell.type"), c("cell.type")),
#'                  cell.type.variables=list(list("iPSC"), list("Endoderm")),
#'                  cell.type.labels=c("iPSC", "Endoderm"),
#'                  n.cells=2,
#'                  features=tran_ids,
#'                  point.size=0.5,
#'                  level="splicing",
#'                  event.type=c("SE", "MXE", "RI", "A5SS", "A3SS"),
#'                  seed=1
#'                 )
#'
#' # Check output
#' marvel$PCA$PSI$Results
#' marvel$PCA$PSI$Plot
#' marvel$PCA$PSI$Plot.Elbow

RunPCA <- function(MarvelObject, cell.type.columns, cell.type.variables, cell.type.labels, n.cells, features, point.size=0.5, level, event.type=NULL, seed=NULL,
                   retrieve.non.outliers=FALSE, pc1.min=NULL, pc1.max=NULL, pc2.min=NULL, pc2.max=NULL,
                   remove.outliers=FALSE
                   ) {

    
    if(level=="splicing") {
        
        RunPCA.PSI(MarvelObject=MarvelObject,
                   cell.type.columns=cell.type.columns,
                   cell.type.variables=cell.type.variables,
                   cell.type.labels=cell.type.labels,
                   n.cells=n.cells,
                   features=features,
                   point.size=point.size,
                   event.type=event.type,
                   seed=seed,
                   retrieve.non.outliers=retrieve.non.outliers,
                   pc1.min=pc1.min,
                   pc1.max=pc1.max,
                   pc2.min=pc2.min,
                   pc2.max=pc2.max,
                   remove.outliers=remove.outliers
                   )

    } else if(level=="gene") {
        
        RunPCA.Exp(MarvelObject=MarvelObject,
                   cell.type.columns=cell.type.columns,
                   cell.type.variables=cell.type.variables,
                   cell.type.labels=cell.type.labels,
                   n.cells=n.cells,
                   features=features,
                   point.size=point.size,
                   retrieve.non.outliers=retrieve.non.outliers,
                   pc1.min=pc1.min,
                   pc1.max=pc1.max,
                   pc2.min=pc2.min,
                   pc2.max=pc2.max,
                   remove.outliers=remove.outliers
                   )
        
    }
    
}
