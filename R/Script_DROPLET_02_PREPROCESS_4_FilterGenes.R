#' @title Filter specific gene types
#'
#' @description Retain genes of specific type, e.g., protein-coding genes.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{AnnotateGenes.10x} function.
#' @param gene.type Character string. Gene type to keep. Specification should match that of GTF.
#'
#' @return An object of class S3 containing the updated slots \code{MarvelObject$gene.metadata}, \code{MarvelObject$gene.norm.matrix}, \code{MarvelObject$sj.metadata}, and \code{MarvelObject$sj.count.matrix}.
#'
#' @importFrom plyr join
#' @import Matrix
#'
#' @export
#'
#' @examples
#'
#' # Load un-processed MARVEL object
#' marvel.demo.10x.raw <- readRDS(system.file("extdata/data",
#'                                "marvel.demo.10x.raw.rds",
#'                                package="MARVEL")
#'                                )
#'
#' # Annotate gene metadata
#' marvel.demo.10x <- AnnotateGenes.10x(MarvelObject=marvel.demo.10x.raw)
#'
#' # Annotate junction metadata
#' marvel.demo.10x <- AnnotateSJ.10x(MarvelObject=marvel.demo.10x)
#'
#' # Validate junctions
#' marvel.demo.10x <- ValidateSJ.10x(MarvelObject=marvel.demo.10x)
#'
#' # Subset CDS genes
#' marvel.demo.10x <- FilterGenes.10x(MarvelObject=marvel.demo.10x,
#'                           gene.type="protein_coding"
#'                           )

FilterGenes.10x <- function(MarvelObject, gene.type="protein_coding") {
    
    # Define arguments
    MarvelObject <- MarvelObject
    gene.metadata <- MarvelObject$gene.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    sj.metadata <- MarvelObject$sj.metadata
    df.sj.count <- MarvelObject$sj.count.matrix
    gene_type <- gene.type
    
    # Example arguments
    #MarvelObject <- marvel
    #gene.metadata <- MarvelObject$gene.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #sj.metadata <- MarvelObject$sj.metadata
    #df.sj.count <- MarvelObject$sj.count.matrix
    #gene_type <- c("protein_coding")
    
    #########################################################
    
    # Report number of genes to be retained
    index <- which(gene.metadata$gene_type %in% gene_type)
    
    message(paste(length(index), " of ", nrow(gene.metadata),  " genes met filtering criteria and retained", sep=""))
    
    # Subset gene metadata
    gene.metadata <- gene.metadata[index, ]
    
    # Subset normalised gene matrix
    df.gene.norm <- df.gene.norm[gene.metadata$gene_short_name, ]
    table(rownames(df.gene.norm)==gene.metadata$gene_short_name)
    
    # Report number of SJ to be retained
    index <- which(sj.metadata$gene_short_name.start %in% gene.metadata$gene_short_name)
    
    message(paste(length(index), " of ", nrow(sj.metadata),  " SJ met filtering criteria and retained", sep=""))
    
    # Subset gene metadata
    sj.metadata <- sj.metadata[index, ]
    
    # Subset normalised gene matrix
    df.sj.count <- df.sj.count[sj.metadata$coord.intron, ]
    table(rownames(df.sj.count)==sj.metadata$coord.intron)
    
    #########################################################
    
    # Update slots
    MarvelObject$gene.metadata <- gene.metadata
    MarvelObject$gene.norm.matrix <- df.gene.norm
    MarvelObject$sj.metadata <- sj.metadata
    MarvelObject$sj.count.matrix <- df.sj.count
 
    message("Normalised gene and SJ metadata and matrix updated")
 
    # Return final object
    return(MarvelObject)
        
}
