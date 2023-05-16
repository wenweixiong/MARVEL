#' @title Preprocess rMATS skipped-exon (SE) coordinates
#'
#' @description Preprocess rMATS skipped-exon (SE) coordinates for MARVEL input and annotate the gene type for each gene using the GTF provided.
#'
#' @param file Data frame. Tab-delimited file output from rMATS. Typically the file is named \code{fromGTF.SE.txt}.
#' @param GTF Data frame. The same GTF file used in the rMATS step.
#'
#' @return An object of class data frame that may be used as input for MARVEL
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export

Preprocess_rMATS.SE <- function(file, GTF) {
    
    # Define arguments
    df <- file
    
    # Example arguments
    #path <- "/Users/seanwen/Documents/U2AF1_2019/Linker/rMATS/ASEvents/"
    #file <- "fromGTF.SE.txt"
    #df <- read.table(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    #path <- "/Users/seanwen/Documents/U2AF1_2019/GTF/"
    #file <- "gencode.v31.annotation.gtf"
    #gtf <- as.data.frame(data.table::fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE, quote=""))
    
    ###############################################################
    
    # Remove "X" prefix from column names
    names(df) <- gsub("^X", "", names(df))
    
    # Create tran_id: +ve strand
        # Subset
        . <- df[which(df$strand=="+"), ]
        
        # Convert coordinates to reflect exon
        .$exonStart_0base <- .$exonStart_0base + 1
        .$upstreamES <- .$upstreamES + 1
        .$downstreamES <- .$downstreamES + 1
        
        # Create tran_od
        .$tran_id <- paste(.$chr, ":", .$upstreamES, ":", .$upstreamEE, ":+@",
                           .$chr, ":", .$exonStart_0base, ":", .$exonEnd, ":+@",
                           .$chr, ":", .$downstreamES, ":", .$downstreamEE,
                           sep=""
                           )
                           
        # Save as new object
        df.pos <- .
                           
    # Create tran_id: -ve strand
        # Subset
        . <- df[which(df$strand=="-"), ]
        
        # Convert coordinates to reflect exon
        .$exonStart_0base <- .$exonStart_0base + 1
        .$upstreamES <- .$upstreamES + 1
        .$downstreamES <- .$downstreamES + 1
        
        # Create tran_od
        .$tran_id <- paste(.$chr, ":", .$downstreamES, ":", .$downstreamEE, ":-@",
                           .$chr, ":", .$exonStart_0base, ":", .$exonEnd, ":-@",
                           .$chr, ":", .$upstreamES, ":", .$upstreamEE,
                           sep=""
                           )
                           
        # Save as new object
        df.neg <- .
                           
    # Merge
    df <- rbind.data.frame(df.pos, df.neg)

    # Subset relavent columns
    df <- df[, c("tran_id", "GeneID", "geneSymbol")]
    names(df)[c(2:3)] <- c("gene_id", "gene_short_name")

    # Keep unique entries
    df <- unique(df)

    # Annotate with gene_type
        # Build gene reference table
            # Subset genes
            ref <- gtf[which(gtf$V3=="gene"), ]
            
            # Subset selected attributes
                # gene_id
                gene_id <- strsplit(ref$V9, split=";")
                gene_id <- sapply(gene_id, function(x) grep("gene_id", x, value=TRUE))
                gene_id <- gsub("gene_id", "", gene_id)
                gene_id <- gsub(" ", "", gene_id)
                gene_id <- gsub("\"", "", gene_id)
                head(gene_id)
                
                # gene_type
                gene_type <- strsplit(ref$V9, split=";")
                gene_type <- sapply(gene_type, function(x) grep("gene_type", x, value=TRUE))
                gene_type <- gsub("gene_type", "", gene_type)
                gene_type <- gsub(" ", "", gene_type)
                gene_type <- gsub("\"", "", gene_type)
                head(gene_type)

                # Create new columns
                ref$gene_id <- gene_id
                ref$gene_type <- gene_type

                # Keep unique entries
                ref <- unique(ref[, c("gene_id", "gene_type")])
                
        # Annotate with attributes
        df <- join(df, ref, by="gene_id", type="left")

    # Collapse duplicate entries
        # Tabulate freq
        . <- as.data.frame(table(df$tran_id))
        names(.) <- c("tran_id", "freq")
        tran_id.unique <- as.character(.[which(.$freq == 1), "tran_id"])
        tran_id.dup <- as.character(.[which(.$freq > 1), "tran_id"])
        
        if(length(tran_id.dup) != 0) {
        
            # Split data frame
            df.unique <- df[which(df$tran_id %in% tran_id.unique), ]
            df.dup <- df[which(df$tran_id %in% tran_id.dup), ]
            
            # Collapse duplicates
            tran_ids <- unique(df.dup$tran_id)
            
            .list <- list()
            
            for(i in 1:length(tran_ids)) {
            
                . <- df.dup[which(df.dup$tran_id == tran_ids[i]), ]
                .list[[i]] <- data.frame("tran_id"=tran_ids[i],
                                         "gene_id"=paste(.$gene_id, collapse="|"),
                                         "gene_short_name"=paste(.$gene_short_name, collapse="|"),
                                         "gene_type"=paste(.$gene_type, collapse="|"),
                                         stringsAsFactors=FALSE
                                         )
                                         
            }
            
            df.dup <- do.call(rbind.data.frame, .list)
            
            # Merge
            df <- rbind.data.frame(df.unique, df.dup)
            
        }
        
    # Check if gene types found for all genes
    if(sum(is.na(df$gene_type) != 0)) {
        
        "Not all genes were successfully annotated with gene_type from GTF file. Please check that the GTF version/file provided is the same as that used in your rMATS step. This is a warning message, not an error."
        
        
    }
   
   ###############################################################
   
   return(df)

}
