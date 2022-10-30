#' @title Pre-flight check
#'
#' @description Ensures only overlapping cells found in both gene and splice junction data are retained. Also ensures matrix columns matches cell IDs in sample metadata and matrix rows matches gene name or splice junction coordinates in feature metadata.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{FilterGenes.10x} function.
#'
#' @return An object of class S3 containing updated slots \code{MarvelObject$gene.norm.matrix}, \code{MarvelObject$sample.metadata}, \code{MarvelObject$gene.metadata}, \code{MarvelObject$gene.count.matrix}, \code{MarvelObject$sj.count.matrix}, \code{MarvelObject$sj.metadata}.
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
#'
#' # Pre-flight check
#' marvel.demo.10x <- CheckAlignment.10x(MarvelObject=marvel.demo.10x)

CheckAlignment.10x <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df.gene.norm <- MarvelObject$gene.norm.matrix
    sample.metadata <- MarvelObject$sample.metadata
    gene.metadata <- MarvelObject$gene.metadata
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    sj.metadata <- MarvelObject$sj.metadata
    
    # Example arguments
    #MarvelObject <- marvel
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #sample.metadata <- MarvelObject$sample.metadata
    #gene.metadata <- MarvelObject$gene.metadata
    #df.gene.count <- MarvelObject$gene.count.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #sj.metadata <- MarvelObject$sj.metadata

    #########################################################
    
    # Subset overlapping genes in normalised matrix and raw count matrix
        # Report progress
        message("Matching gene list in normalised and raw (count) gene matrices...")
        
        # Find overlaps
        overlap <- intersect(rownames(df.gene.norm), row.names(df.gene.count))
    
        # Report numbers
        message(paste(length(rownames(df.gene.norm)), " genes found in normalised gene matrix", sep=""))
        message(paste(length(rownames(df.gene.count)), " genes found in raw (count) gene matrix", sep=""))
        message(paste(length(overlap), " overlapping genes found and retained", sep=""))
        
        # Subset normalised matrix
        df.gene.norm <- df.gene.norm[overlap, ]
        
        # Subset raw (count) matrix
        df.gene.count <- df.gene.count[overlap, ]
        
        # Subset gene metadata
        gene.metadata <- gene.metadata[which(gene.metadata$gene_short_name %in% overlap), ]
   
   # Subset genes in which at least 1 SJ is found
        # Find overlaps
        overlap <- intersect(gene.metadata$gene_short_name, sj.metadata$gene_short_name.start)
      
        # Report numbers
        paste(length(overlap), " genes with at least 1 SJ found and retained", sep="")
        
        # Subset normalised gene matrix
        df.gene.norm <- df.gene.norm[overlap, ]
        
        # Subset raw (count) gene matrix
        df.gene.count <- df.gene.count[overlap, ]
        
        # Subset gene metadata
        gene.metadata <- gene.metadata[which(gene.metadata$gene_short_name %in% overlap), ]
        
        # Subset SJ metadata
        sj.metadata <- sj.metadata[which(sj.metadata$gene_short_name.start %in% overlap), ]
        
        # Subset raw (count) gene matrix
        df.sj.count <- df.sj.count[sj.metadata$coord.intron, ]
        
    # Subset overlapping cells across all 3 matrices
        # Report progress
        message("Matching cells across normalised gene, raw (count) gene, raw (count) SJ matrices...")
        
        # Find overlaps
        .list <- list(colnames(df.gene.norm),
                      colnames(df.gene.count),
                      colnames(df.sj.count)
                      )
        
        overlap <- Reduce(intersect, .list)
        
        # Report numbers
        message(paste(length(colnames(df.gene.norm)), " cells found in normalised gene matrix", sep=""))
        message(paste(length(colnames(df.gene.count)), " cells found in raw (count) gene matrix", sep=""))
        message(paste(length(colnames(df.sj.count)), " cells found in raw (count) SJ matrix", sep=""))
        message(paste(length(overlap), " overlapping cells found and retained", sep=""))
        
        # Subset normalised gene matrix
        df.gene.norm <- df.gene.norm[, overlap]
        
        # Subset raw (count) gene matrix
        df.gene.count <- df.gene.count[, overlap]
        
        # Subset raw (count) SJ matrix
        df.sj.count <- df.sj.count[, overlap]
        
        # Subset sample metadata
        sample.metadata <- sample.metadata[which(sample.metadata$cell.id %in% overlap), ]
        
    # Check column alignment
        # Report progress
        message("Checking column (sample) alignment...")
        
        # Sample metadata vs normalised gene expression matrix
        index.l <- table(sample.metadata$cell.id==colnames(df.gene.norm))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message("Sample metadata and normalised gene matrix column names MATCHED")
            
        } else {
            
            
            message("Sample metadata and normalised gene matrix column names NOT MATCHED")

        }
        
        # Normalised vs raw (count) gene expression matrix
        index.l <- table(colnames(df.gene.norm)==colnames(df.gene.count))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message("Normalised and raw (count) gene matrix column names MATCHED")
            
        } else {
            
            
            message("Normalised and raw (count) gene matrix column names NOT MATCHED")

        }
        
        # Raw (count) gene and SJ expression matrix
        index.l <- table(colnames(df.gene.count)==colnames(df.sj.count))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message("Raw (Count) gene and SJ matrix column names MATCHED")
            
        } else {
            
            
            message("Raw (Count) gene and SJ matrix column names NOT MATCHED")

        }

    # Check row alignment
        # Report progress
        message("Checking row (gene names/SJ coordinates) alignment...")
        
        # Gene metadata vs normalised gene expression matrix
        index.l <- table(gene.metadata$gene_short_name==rownames(df.gene.norm))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message("Gene metadata and normalised gene matrix row names MATCHED")
            
        } else {
            
            
            message("Gene metadata and normalised gene matrix row names NOT MATCHED")

        }
        
        # Normalised vs raw (count) gene expression matrix
        index.l <- table(rownames(df.gene.norm)==rownames(df.gene.count))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message("Normalised and raw (count) gene matrix row names MATCHED")
            
        } else {
            
            
            message("Normalised and raw (count) gene matrix row names NOT MATCHED")

        }
    
        # SJ metadata vs raw (count) SJ expression matrix
        index.l <- table(sj.metadata$coord.intron==rownames(df.sj.count))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message("SJ metadata and raw (count) SJ matrix row names MATCHED")
            
        } else {
            
            message("SJ metadata and raw (count) SJ matrix row names NOT MATCHED")

        }
        
    # Report final numbers for downstream analysis
    n.cells <- nrow(sample.metadata)
    n.genes <- nrow(gene.metadata)
    n.sj <- nrow(sj.metadata)
    
    message(paste(n.cells, " cells and ", n.genes, " genes consisting of ", n.sj, " splice junctions included for further analysis", sep=""))
    
    #########################################################
    
    # Update slots
    MarvelObject$gene.norm.matrix <- df.gene.norm
    MarvelObject$sample.metadata <- sample.metadata
    MarvelObject$gene.metadata <- gene.metadata
    MarvelObject$gene.count.matrix <- df.gene.count
    MarvelObject$sj.count.matrix <- df.sj.count
    MarvelObject$sj.metadata <- sj.metadata
        
    # Return final object
    return(MarvelObject)
            
}


