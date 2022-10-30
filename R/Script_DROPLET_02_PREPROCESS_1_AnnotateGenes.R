#' @title Annotate splice junctions
#'
#' @description Annotates the each gene in the gene metadata with the gene type, e.g. protein-coding, antisense etc.. Annotations are retrieved from GTF. Only genes found in gene metadata and GTF will be retained.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject.10x} function.
#'
#' @return An object of class S3 containing the updated slots \code{MarvelObject$gene.metadata} and \code{gene.norm.matrix}.
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

AnnotateGenes.10x <- function(MarvelObject) {
    
    # Define arguments
    MarvelObject <- MarvelObject
    gtf <- MarvelObject$gtf
    gene.metadata <- MarvelObject$gene.metadata
    df.gene.norm <- MarvelObject$gene.norm
    
    # Example arguments
    #MarvelObject <- marvel.demo.10x.raw
    #gtf <- MarvelObject$gtf
    #gene.metadata <- MarvelObject$gene.metadata
    #df.gene.norm <- MarvelObject$gene.norm
    
    ###############################################################
    
    # Subset gene entries
    gtf <- gtf[which(gtf$V3=="gene"), ]

    # Retrieve selected attributes
        # gene names
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("gene_name", x, value=TRUE))
        . <- gsub("gene_name", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)

        gtf$gene_short_name <- .

        # gene type
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("gene_biotype", x, value=TRUE))
        . <- gsub("gene_biotype", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)

        gtf$gene_type <- .
    
        # Subset unique entries
        gtf <- gtf[,c("gene_short_name", "gene_type")]
        gtf <- unique(gtf)
        
        # Remove genes with >1 annotation
        #. <- as.data.frame(table(gtf$gene_short_name))
        
    # Subset overlapping genes
    overlap <- intersect(gene.metadata$gene_short_name, unique(gtf$gene_short_name))
    message(paste(nrow(gene.metadata), " genes found in gene metadata", sep=""))
    message(paste(length(overlap), " genes overlapped with GTF and will be subset-ed", sep=""))
    
    gene.metadata <- gene.metadata[which(gene.metadata$gene_short_name %in% overlap), , drop=FALSE]
    df.gene.norm <- df.gene.norm[gene.metadata$gene_short_name, ]
    
    # Annotate gene metadata
    gene.metadata <- join(gene.metadata, gtf, by="gene_short_name", type="left", match="first")
    sum(is.na(gene.metadata$gene_type))
    
    message("Gene type annotated")
    
    # Check alignment
    table(rownames(df.gene.norm)==gene.metadata$gene_short_name)
    
    ###############################################################
    
    # Update slots
    MarvelObject$gene.metadata <- gene.metadata
    MarvelObject$gene.norm.matrix <- df.gene.norm
    
    message("Gene metadata and normalised gene expression matrix updated")
        
    # Return final object
    return(MarvelObject)
        
}
