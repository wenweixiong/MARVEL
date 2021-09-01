#' @title Parse Gene Transfer File (GTF)
#'
#' @description
#' \code{ParseGTF} parses the gene transfer file (GTF) for downstream nonsense-mediated decay (NMD) prediction.
#'
#' @details
#' This function parses the GTF in order to generate new columns for gene IDs, transcript IDs, and transcript type. These information are extracted from the attribute (9th) column for a standard GTF. These information will be used for downstream NMD prediction.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues.PSI} function.
#' @param gtf Data frame. Nine-column data frame containing gene, transcript, and exon information. Please see https://www.gencodegenes.org/pages/data_format.html for formatting. This GTF file should be the same as the one used for detecting the initial alternative splicing events, e.g. rMATS etc.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$NMD$GTF}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @import textclean
#'
#' @examples
#' # Load input
#' path_to_file <- system.file("extdata/Data", "gencode.v31.annotation.gtf", package="MARVEL")
#' gtf <- read.table(path_to_file, sep="\t", header=FALSE, stringsAsFactors=FALSE,
#'              na.strings="NA")
#'
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- ParseGTF(MarvelObject=marvel, gtf=gtf)
#'
#' # Check output
#' marvel$NMD$GTF[1:5, ]

ParseGTF <- function(MarvelObject, gtf) {

    # Define arguments
    gtf <- gtf

    # Parse attributes
    print("Parsing attribute column...")
    
    attr <- strsplit(gtf$V9, split=";")

    # Retrieve selected attributes
        # gene_id
        print("Retrieving gene_id...")
        
        . <- sapply(attr, function(x) grep("gene_id", x, value=TRUE))
        gtf$gene_id <- mgsub(., c("gene_id", " ", "\""), "")
        
        # transcript_id
        print("Retrieving transcript_id...")
        
        . <- sapply(attr, function(x) grep("transcript_id", x, value=TRUE))
        gtf$transcript_id <- mgsub(., c("transcript_id", " ", "\""), "")
        
        # transcript_type
        print("Retrieving transcript_type...")
        
        . <- sapply(attr, function(x) grep("transcript_type", x, value=TRUE))
        gtf$transcript_type <- mgsub(., c("transcript_type", " ", "\""), "")
        
    # Remove attribute column
    gtf$V9 <- NULL

    # Save to new slot
    MarvelObject$NMD$GTF <- gtf
    
    # Return MARVEL object
    return(MarvelObject)

}
