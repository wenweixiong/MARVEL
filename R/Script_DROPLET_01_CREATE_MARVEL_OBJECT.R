#' @title Create Marvel object for droplet-based RNA-sequencing data
#'
#' @description Creates an S3 object named \code{Marvel} for downstream analysis, specifically for droplet-based RNA-sequencing data.
#'
#' @param gene.norm.matrix Sparse matrix. UMI-collapsed, normalised, non-log2-transformed gene expression matrix.
#' @param gene.norm.pheno Data frame. Sample metadata for annotating \code{gene.norm.matrix} columns with cell IDs.
#' @param gene.norm.feature Data frame. Gene metadata for annotating \code{gene.norm.matrix} rows with gene names.
#' @param gene.count.matrix Sparse matrix. UMI-collapsed, non-normalised (raw counts), non-log2-transformed gene expression matrix.
#' @param gene.count.pheno Data frame. Sample metadata for annotating \code{gene.count.matrix} columsn with cell IDs.
#' @param gene.count.feature Data frame. Gene metadata for annotating \code{gene.count.matrix} rows with gene names.
#' @param sj.count.matrix Sparse matrix. UMI-collapsed, non-normalised (raw counts), non-log2-transformed splice junction expression matrix.
#' @param sj.count.pheno Data frame. Sample metadata for annotating \code{sj.count.matrix} columsn with cell IDs.
#' @param sj.count.feature Data frame. Splice junction metadata for annotating \code{sj.count.matrix} rows with splice junction coordinates.
#' @param pca Data frame. Coordinates of PCA/tSNE/UMAP.
#' @param gtf Data frame. GTF used in cellranger. Will be used for annotating splice junctions downstream.
#'
#' @return An object of class S3.
#'
#' @importFrom plyr join
#' @import Matrix
#'
#' @export
#'
#' @examples
#' # Retrieve, observe format of pre-saved input files
#' marvel.demo.10x.raw <- readRDS(system.file("extdata/data",
#'                                "marvel.demo.10x.raw.rds",
#'                                package="MARVEL")
#'                                )
#' # Gene expression (Normalised)
#'     # Matrix
#'     df.gene.norm <- marvel.demo.10x.raw$gene.norm.matrix
#'     df.gene.norm[1:5, 1:5]
#'
#'     # phenoData
#'     df.gene.norm.pheno <- marvel.demo.10x.raw$sample.metadata
#'     head(df.gene.norm.pheno)
#'
#'     # featureData
#'     df.gene.norm.feature <- data.frame("gene_short_name"=rownames(df.gene.norm),
#'                                        stringsAsFactors=FALSE
#'                                        )
#'     head(df.gene.norm.feature)
#'
#' # Gene expression (Counts)
#'     # Matrix
#'     df.gene.count <- marvel.demo.10x.raw$gene.count.matrix
#'     df.gene.count[1:5, 1:5]
#'
#'     # phenoData
#'     df.gene.count.pheno <- data.frame("cell.id"=colnames(df.gene.count),
#'                                        stringsAsFactors=FALSE
#'                                        )
#'     head(df.gene.count.pheno)
#'
#'     # featureData
#'     df.gene.count.feature <- data.frame("gene_short_name"=rownames(df.gene.count),
#'                                        stringsAsFactors=FALSE
#'                                        )
#'     head(df.gene.count.feature)
#'
#' # SJ (Counts)
#'     # Matrix
#'     df.sj.count <- marvel.demo.10x.raw$sj.count.matrix
#'     df.sj.count[1:5, 1:5]
#'
#'     # phenoData
#'     df.sj.count.pheno <- data.frame("cell.id"=colnames(df.sj.count),
#'                                      stringsAsFactors=FALSE
#'                                      )
#'     head(df.sj.count.pheno)
#'
#'     # featureData
#'     df.sj.count.feature <- data.frame("coord.intron"=rownames(df.sj.count),
#'                                        stringsAsFactors=FALSE
#'                                        )
#'     head(df.sj.count.feature)
#'
#' # tSNE coordinates
#' df.coord <- marvel.demo.10x.raw$pca
#' head(df.coord)
#'
#' # GTF
#' gtf <- marvel.demo.10x.raw$gtf
#' head(gtf)
#'
#' # Create MARVEL object
#' marvel.demo.10x <- CreateMarvelObject.10x(gene.norm.matrix=df.gene.norm,
#'                      gene.norm.pheno=df.gene.norm.pheno,
#'                      gene.norm.feature=df.gene.norm.feature,
#'                      gene.count.matrix=df.gene.count,
#'                      gene.count.pheno=df.gene.count.pheno,
#'                      gene.count.feature=df.gene.count.feature,
#'                      sj.count.matrix=df.sj.count,
#'                      sj.count.pheno=df.sj.count.pheno,
#'                      sj.count.feature=df.sj.count.feature,
#'                      pca=df.coord,
#'                      gtf=gtf
#'                      )

CreateMarvelObject.10x <- function(gene.norm.matrix=NULL,
                                   gene.norm.pheno=NULL,
                                   gene.norm.feature=NULL,
                                   gene.count.matrix=NULL,
                                   gene.count.pheno=NULL,
                                   gene.count.feature=NULL,
                                   sj.count.matrix=NULL,
                                   sj.count.pheno=NULL,
                                   sj.count.feature=NULL,
                                   pca=NULL,
                                   gtf=NULL
                                   ) {
        
    # Annotate row and columns for matrices
        # Track progress
        message("Annotating rows and columns for matrices...")
        
        # Normalised gene expression matrix
        colnames(gene.norm.matrix) <- gene.norm.pheno$cell.id
        rownames(gene.norm.matrix) <- gene.norm.feature$gene_short_name
        
        # Raw (counts) gene expression matrix
        colnames(gene.count.matrix) <- gene.count.pheno$cell.id
        rownames(gene.count.matrix) <- gene.count.feature$gene_short_name
        
        # Raw (counts) SJ matrix
        colnames(sj.count.matrix) <- sj.count.pheno$cell.id
        rownames(sj.count.matrix) <- sj.count.feature$coord.intron
        
    # Create s3 object
    s3 <- list()
    class(s3) <- "Marvel"
    
    # Fill slots
    message("Filling in slots for MARVEL object...")
    
    s3$sample.metadata <- gene.norm.pheno
    s3$gene.metadata <- gene.norm.feature
    s3$gene.norm.matrix <- gene.norm.matrix
    s3$gene.count.matrix <- gene.count.matrix
    s3$sj.count.matrix <- sj.count.matrix
    s3$pca <- pca
    s3$gtf <- gtf
    
    # Returen final object
    message("MARVEL object created")
    return(s3)
        
}


