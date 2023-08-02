#' @title Preprocess rMATS retained intron (RI) coordinates for Bedtools
#'
#' @description Preprocess rMATS retained intron (RI) coordinates into a BED file suitable for counting intronic coverage with Bedtools
#'
#' @param file Data frame. Tab-delimited file output from rMATS. Typically the file is named \code{fromGTF.RI.txt}.
#'
#' @return An object of class data frame that may be used as input for MARVEL
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export

PrepareBedFile.RI <- function(file) {
    
    # Define arguments
    df <- file
    
    # Example arguments
    path <- "/Users/seanwen/Documents/U2AF1_2019/Linker/rMATS/ASEvents/"
    file <- "fromGTF.RI.txt"
    df <- read.table(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)

    ###############################################################
    
    # Remove "X" prefix from column names
    names(df) <- gsub("^X", "", names(df))
    
    # Subset intron coordinates: +ve strand
        # Subset
        . <- df[which(df$strand=="+"), ]
        
        # Convert coordinates to reflect exon
        #.$riExonStart_0base <- .$riExonStart_0base + 1
        #.$downstreamES <- .$downstreamES + 1

        # Subset coordinates
        . <- .[,c("chr", "upstreamEE", "downstreamES")]
        
        # Save as new object
        df.pos <- .
                           
    # Create tran_id: -ve strand
        # Subset
        . <- df[which(df$strand=="-"), ]
        
        # Convert coordinates to reflect exon
        #.$riExonStart_0base <- .$riExonStart_0base + 1
        #.$downstreamES <- .$downstreamES + 1

        # Subset coordinates
        . <- .[,c("chr", "upstreamEE", "downstreamES")]
        
        # Save as new object
        df.neg <- .
                           
    # Merge
    df <- rbind.data.frame(df.pos, df.neg)

    # Remove exon coordinates
    #df$upstreamEE <- df$upstreamEE + 1
    #df$downstreamES <- df$downstreamES - 1

    # Reorder
    table(df$chr)
    df$chr <- factor(df$chr,
                    levels=c("chr1", "chr2", "chr3", "chr4", "chr5",
                             "chr6", "chr7", "chr8", "chr9", "chr10",
                             "chr11", "chr12", "chr13", "chr14", "chr15",
                             "chr16", "chr17", "chr18", "chr19", "chr20",
                             "chr21", "chr22", "chrX"
                              ))
    df <- df[order(df$chr, df$upstreamEE),]

    # Check coordinates
    table(df$downstreamES > df$upstreamEE)

    # Keep unique entries
    df <- unique(df)
   
   ###############################################################
   
   return(df)

}
