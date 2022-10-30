#' @title Parse gene transfer file (GTF)
#'
#' @description Parses the gene transfer file (GTF) for downstream nonsense-mediated decay (NMD) prediction.
#'
#' @details
#' This function parses the GTF in order to generate new columns for gene IDs, transcript IDs, and transcript type. These information are extracted from the attribute (9th) column for a standard GTF. These information will be used for downstream NMD prediction.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.PSI} function.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$NMD$GTF}.
#'
#' @importFrom plyr join
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- ParseGTF(MarvelObject=marvel.demo)

ParseGTF <- function(MarvelObject) {

    # Define arguments
    gtf <- MarvelObject$GTF

    # Parse attributes
    message("Parsing attribute column...")
    
    attr <- strsplit(gtf$V9, split=";")

    # Retrieve selected attributes
        # gene_id
        message("Retrieving gene_id...")
        
        . <- sapply(attr, function(x) grep("gene_id", x, value=TRUE))
        gtf$gene_id <- textclean::mgsub(., c("gene_id", " ", "\""), "")
        
        # transcript_id
        message("Retrieving transcript_id...")
        
        . <- sapply(attr, function(x) grep("transcript_id", x, value=TRUE))
        gtf$transcript_id <- textclean::mgsub(., c("transcript_id", " ", "\""), "")
        
        # transcript_type
        message("Retrieving transcript_type...")
        
        . <- sapply(attr, function(x) grep("transcript_type", x, value=TRUE))
        gtf$transcript_type <- textclean::mgsub(., c("transcript_type", " ", "\""), "")
        
    # Remove attribute column
    gtf$V9 <- NULL

    # Save to new slot
    MarvelObject$NMD$GTF <- gtf
    
    # Return MARVEL object
    return(MarvelObject)

}
