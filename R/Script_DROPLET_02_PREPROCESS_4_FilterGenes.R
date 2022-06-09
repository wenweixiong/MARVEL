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
#'
#' @export

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
    
    print(paste(length(index), " of ", nrow(gene.metadata),  " genes met filtering criteria and retained", sep=""))
    
    # Subset gene metadata
    gene.metadata <- gene.metadata[index, ]
    
    # Subset normalised gene matrix
    df.gene.norm <- df.gene.norm[gene.metadata$gene_short_name, ]
    table(rownames(df.gene.norm)==gene.metadata$gene_short_name)
    
    # Report number of SJ to be retained
    index <- which(sj.metadata$gene_short_name.start %in% gene.metadata$gene_short_name)
    
    print(paste(length(index), " of ", nrow(sj.metadata),  " SJ met filtering criteria and retained", sep=""))
    
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
 
    print("Normalised gene and SJ metadata and matrix updated")
 
    # Return final object
    return(MarvelObject)
        
}
