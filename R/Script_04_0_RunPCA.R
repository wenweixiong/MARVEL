#' @title Principle Component Analysis
#'
#' @description
#' \code{RunPCA} performs principle component analysis on splicing and gene expression data.
#'
#' @details
#' This function performs principle component analysis on splicing and gene expression data and visualise cells on a reducted dimension space, i.e. 2D scatterplot.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.types Character string. To indicate which groups of cells that will be used for analysis. Group names should match those in \code{cell.type} column of \code{$SplicePheno} or \code{GenePheno} slot for splicing or gene expression data, respectively.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event or gene for the event or gene, respectively, to be included for analysis.
#' @param features Character string. Vector of \code{tran_id}s or \code{gene_id}s for analysis. Should match \code{tran_id} or \code{gene_id} column of \code{$ValidatedSpliceFeature} or \code{GeneFeature} slot when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param level Character string. Indicate \code{"splicing"} or \code{"gene"} for splicing or gene expression analysis, respectively
#' @param event.type Character string. Only applicable when \code{level} set to \code{"splicing"}.  Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param seed Numeric value. Only applicable when \code{level} set to \code{"splicing"}.  Ensures imputed values for NA PSIs are reproducible.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$PCA$PSI} or \code{MarvelObject$PCA$Gene} when \code{level} set to \code{"splicing"} or \code{"gene"}, respectively. Contains 2D scatterplot in \code{MarvelObject$PCA$PSI$Plot} or \code{MarvelObject$PCA$Gene$Plot} and the corresponding x- and y-coordinates for each sample in \code{MarvelObject$PCA$PSI$Results} or \code{MarvelObject$PCA$Gene$Results}, respectively.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import FactoMineR
#' @import factoextra
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' features <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
#' features <- features$tran_id
#'
#' marvel <- RunPCA(MarvelObject=marvel,
#'                  cell.types="all",
#'                  n.cells=3,
#'                  features=features,
#'                  point.size=2.5,
#'                  level="splicing",
#'                  event.type="all",
#'                  seed=1
#'                  )
#'
#' marvel$PCA$PSI$Results
#' marvel$PCA$PSI$Plot

RunPCA <- function(MarvelObject, cell.types, n.cells, features, point.size, level, event.type=NULL, seed=NULL) {

    
    if(level=="splicing") {
        
        RunPCA.PSI(MarvelObject=MarvelObject,
                   cell.types=cell.types,
                   n.cells=n.cells,
                   features=features,
                   point.size=point.size,
                   event.type=event.type,
                   seed=seed
                   )

    } else if(level=="gene") {
        
        RunPCA.Exp(MarvelObject=MarvelObject,
                    cell.types=cell.types,
                    n.cells=n.cells,
                    features=features,
                    point.size=point.size
                   )
        
    }
    
}
