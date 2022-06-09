#' @title Validate splice junctions
#'
#' @description Retains splice junctions whose start and end belong to the same gene.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{AnnotateSJ.10x} function.
#'
#' @return An object of class S3 containing the updated slots \code{MarvelObject$sj.metadata} and \code{MarvelObject$sj.count.matrix}.
#'
#' @importFrom plyr join
#'
#' @export

ValidateSJ.10x <- function(MarvelObject) {
    
    # Define arguments
    MarvelObject <- MarvelObject
    sj.metadata<- MarvelObject$sj.metadata
    df.sj.count <- MarvelObject$sj.count.matrix
    
    # Example arguments
    #MarvelObject <- marvel
    #sj.metadata<- MarvelObject$sj.metadata
    #df.sj.count <- MarvelObject$sj.count.matrix
    
    #########################################################
            
    # Report SJ type
        # Annotated + unambiguous
        sj.types.1 <- c("start_known.single.gene|end_known.single.gene|same")
        
        index <- which(sj.metadata$sj.type %in% sj.types.1)
        
        if(length(index) >= 1) {
            
            print(paste(length(index), " annotated and uniquely-mapped SJ identified and retained", sep=""))
            
        }
        
        # Unannotated
        sj.types.2 <- c("start_known.single.gene|end_unknown.gene",
                        "start_unknown.gene|end_known.single.gene",
                        "start_unknown.gene|end_unknown.gene"
                        )
                        
        index <- which(sj.metadata$sj.type %in% sj.types.2)
                        
        if(length(index) >= 1) {
            
            print(paste(length(index), " unannotated SJ identified and removed", sep=""))
            
        }
        
        # Ambiguous
        sj.types.3 <- c("start_known.multi.gene|end_known.multi.gene",
                        "start_known.multi.gene|start_known.single.gene",
                        "start_known.single.gene|end_known.multi.gene",
                        "start_known.single.gene|end_known.single.gene|different"
                        )
        
        index <- which(sj.metadata$sj.type %in% sj.types.3)
                        
        if(length(index) >= 1) {
            
            print(paste(length(index), " multi-mapping SJ identified and removed", sep=""))
            
        }
        
        # Unannotated + ambiguous
        sj.types.4 <- c("start_known.multi.gene|end_unknown.gene",
                        "start_unknown.gene|end_known.multi.gene"
                        )
        
        index <- which(sj.metadata$sj.type %in% sj.types.4)
                        
        if(length(index) >= 1) {
            
            print(paste(length(index), " unannotated and multi-mapping SJ identified and removed", sep=""))
            
        }
        
    # Subset relevant SJ type
    sj.metadata <- sj.metadata[which(sj.metadata$sj.type %in% sj.types.1), , drop=FALSE]
    df.sj.count <- df.sj.count[sj.metadata$coord.intron, ]
    
    # Check alignment with SJ matrix
    table(rownames(df.sj.count)==sj.metadata$coord.intron)
   
    #########################################################
    
    # Update slots
    MarvelObject$sj.metadata <- sj.metadata
    MarvelObject$sj.count.matrix <- df.sj.count
 
    print("SJ metadata and SJ count matrix updated")
 
    # Return final object
    return(MarvelObject)
        
}
