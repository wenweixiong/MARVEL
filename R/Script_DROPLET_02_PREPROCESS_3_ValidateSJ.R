#' @title Validate splice junctions
#'
#' @description Retains splice junctions whose start and end belong to the same gene.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{AnnotateSJ.10x} function.
#' @param keep.novel.sj Logical value. If set to \code{TRUE}, novel splice junctions will be retained for downstream analysis. Novel splice junctions are defined as splice junctions with one end reported in GTF while the other was not reported in GTF. Default value is \code{FALSE}.
#'
#' @return An object of class S3 containing the updated slots \code{MarvelObject$sj.metadata} and \code{MarvelObject$sj.count.matrix}.
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

ValidateSJ.10x <- function(MarvelObject, keep.novel.sj=FALSE) {
    
    # Define arguments
    MarvelObject <- MarvelObject
    sj.metadata <- MarvelObject$sj.metadata
    df.sj.count <- MarvelObject$sj.count.matrix
    keep.novel.sj <- keep.novel.sj
    
    # Example arguments
    #MarvelObject <- marvel
    #sj.metadata<- MarvelObject$sj.metadata
    #df.sj.count <- MarvelObject$sj.count.matrix
    #keep.novel.sj <- TRUE
    
    #########################################################
            
    # Report SJ type
        # Annotated + unambiguous
        sj.types.1 <- c("start_known.single.gene|end_known.single.gene|same")
        
        index <- which(sj.metadata$sj.type %in% sj.types.1)
        
        if(length(index) >= 1) {
            
            message(paste(length(index), " annotated and uniquely-mapped SJ identified", sep=""))
            
        }
        
        # Annotated + unambiguous
        sj.types.2 <- c("start_known.single.gene|end_unknown.gene",
                        "start_unknown.gene|end_known.single.gene"
                        )
                        
        index <- which(sj.metadata$sj.type %in% sj.types.2)
        
        if(length(index) >= 1) {
            
            message(paste(length(index), " annotated (at one end) and uniquely-mapped SJ identified", sep=""))
            
        }
        
        # Unannotated
        sj.types.3 <- c("start_unknown.gene|end_unknown.gene"
                        )
                        
        index <- which(sj.metadata$sj.type %in% sj.types.3)
                        
        if(length(index) >= 1) {
            
            message(paste(length(index), " unannotated SJ identified", sep=""))
            
        }
        
        # Ambiguous
        sj.types.4 <- c("start_known.multi.gene|end_known.multi.gene",
                        "start_known.multi.gene|start_known.single.gene",
                        "start_known.single.gene|end_known.multi.gene",
                        "start_known.single.gene|end_known.single.gene|different"
                        )
        
        index <- which(sj.metadata$sj.type %in% sj.types.4)
                        
        if(length(index) >= 1) {
            
            message(paste(length(index), " multi-mapping SJ identified", sep=""))
            
        }
        
        # Unannotated + ambiguous
        sj.types.5 <- c("start_known.multi.gene|end_unknown.gene",
                        "start_unknown.gene|end_known.multi.gene"
                        )
        
        index <- which(sj.metadata$sj.type %in% sj.types.5)
                        
        if(length(index) >= 1) {
            
            message(paste(length(index), " unannotated and multi-mapping SJ identified", sep=""))
            
        }
        
    # Subset relevant SJ type(s)
    if(keep.novel.sj==TRUE) {
        
        # Subset
        sj.metadata <- sj.metadata[which(sj.metadata$sj.type %in% c(sj.types.1, sj.types.2)), , drop=FALSE]
        
        # Annotate novel end with gene name from known end
        index <- which(sj.metadata$sj.type=="start_known.single.gene|end_unknown.gene")
        sj.metadata$gene_short_name.end[index] <- sj.metadata$gene_short_name.start[index]
        
        index <- which(sj.metadata$sj.type=="start_unknown.gene|end_known.single.gene")
        sj.metadata$gene_short_name.start[index] <- sj.metadata$gene_short_name.end[index]
        
        # Check for missing values
        sum(is.na(sj.metadata$gene_short_name.start))
        sum(is.na(sj.metadata$gene_short_name.end))
        table(sj.metadata$gene_short_name.start==sj.metadata$gene_short_name.end)
        
    } else {
        
        # Subset
        sj.metadata <- sj.metadata[which(sj.metadata$sj.type %in% c(sj.types.1)), , drop=FALSE]
        
    }
    
    df.sj.count <- df.sj.count[sj.metadata$coord.intron, ]
    
    # Check alignment with SJ matrix
    table(rownames(df.sj.count)==sj.metadata$coord.intron)
   
    #########################################################
    
    # Update slots
    MarvelObject$sj.metadata <- sj.metadata
    MarvelObject$sj.count.matrix <- df.sj.count
 
    message("SJ metadata and SJ count matrix updated")
 
    # Return final object
    return(MarvelObject)
        
}
