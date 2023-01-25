#' @title Detect alternative last exons
#'
#' @description Detects alternative last exons from GTF. This is a wrapper function for \code{DetectEvents.ALE.PosStrand} and \code{DetectEvents.ALE.NegStrand} functions.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param min.cells Numeric value. The minimum number of cells in which the gene is expressed for the gene to included for splicing event detected and quantification. To be used in conjunction with \code{min.expr} argument. Default value is \code{50}.
#' @param min.expr Numeric value. The minimum expression value for the gene to be considered to be expressed in a cell. Default value is \code{1}.
#' @param track.progress Logical. If set to \code{TRUE}, progress bar will appear to track the progress of the rate-limiting step of this function, which is the extraction of the final exon-exon junctions. Default value is \code{FALSE}.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$SpliceFeature$ALE}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- DetectEvents.ALE(MarvelObject=marvel.demo,
#'                                 min.cells=5,
#'                                 min.expr=1,
#'                                 track.progress=FALSE
#'                                 )

DetectEvents.ALE <- function(MarvelObject, min.cells=50, min.expr=1, track.progress=FALSE) {

    # Define arguments
    df <- MarvelObject$GTF
    min.cells <- min.cells
    min.expr <- min.expr
    track.progress <- track.progress
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$GTF
    #min.cells <- 27
    #min.expr <- 1
    
    # Check if GTF provided
    if(is.null(df)) {
        
        message("Please provide GTF during creating MARVEL object step to proceed")
        
        return(MarvelObject)
        
    }
    
    # Remove novel transcripts
    index <- grep("MSTRG", df$V9, fixed=TRUE)
    
    if(length(index) != 0) {
        
        df <- df[-index, ]
        
    }
    
    # Parse GTF
    message("Parsing GTF...")
            
    attr <- strsplit(df$V9, split=";")
    . <- sapply(attr, function(x) grep("gene_id", x, value=TRUE))
    df$gene_id <- textclean::mgsub(., c("gene_id", " ", "\""), "")
        
    # Detect, compute PSI on +ve strand
    message("Analysing +ve strand...")
    
    MarvelObject.PosStrand <- DetectEvents.ALE.PosStrand(MarvelObject=MarvelObject,
                              parsed.gtf=df,
                              min.cells=min.cells,
                              min.expr=min.expr,
                              track.progress=track.progress
                              )
        
    # Detect, compute PSI on -ve strand
    message("Analysing -ve strand...")
    
    MarvelObject.NegStrand <- DetectEvents.ALE.NegStrand(MarvelObject=MarvelObject,
                              parsed.gtf=df,
                              min.cells=min.cells,
                              min.expr=min.expr,
                              track.progress=track.progress
                              )
    
    ######################################################################
    ###################### RETURN FINAL OBJECTS ##########################
    ######################################################################

    # Merge +ve, -ve strand
    df.feature.merged <- rbind.data.frame(MarvelObject.PosStrand$SpliceFeature$ALE.PosStrand,
                                          MarvelObject.NegStrand$SpliceFeature$ALE.NegStrand
                                          )

    # Save to new slots
    MarvelObject$SpliceFeature$ALE <- df.feature.merged
    
    # Remove intermediate objects
    MarvelObject.PosStrand <- NULL
    MarvelObject.NegStrand <- NULL
    
    # Remove intermediate slots
    MarvelObject.PosStrand$SpliceFeature$ALE.PosStrand <- NULL
    MarvelObject.NegStrand$SpliceFeature$ALE.NegStrand <- NULL
    
    # Return final object
    return(MarvelObject)

    
}
