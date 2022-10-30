#' @title Annotate splice junctions
#'
#' @description Annotates the splice junctions by assigning the gene name to the start and end of the splice junction. Annotations are retrieved from GTF.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{AnnotateGenes.10x} function.
#'
#' @return An object of class S3 containing the updated slot \code{MarvelObject$sj.metadata}.
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

AnnotateSJ.10x <- function(MarvelObject) {
    
    # Define arguments
    MarvelObject <- MarvelObject
    gtf <- MarvelObject$gtf
    df.sj.count <- MarvelObject$sj.count.matrix
    
    # Example arguments
    #MarvelObject <- marvel
    #gtf <- MarvelObject$gtf
    #df.sj.count <- MarvelObject$sj.count.matrix
    
    #########################################################
    
    # Create SJ metadata
    message("Creating splice junction metadata...")
    
    df <- data.frame("coord.intron"=rownames(df.sj.count), stringsAsFactors=FALSE)
    
    . <- strsplit(rownames(df.sj.count), split=":", fixed=TRUE)
    df$chr <- sapply(., function(x) {as.character(x[1])})
    df$start <- sapply(., function(x) {as.character(x[2])})
    df$end <- sapply(., function(x) {as.character(x[3])})
    
    # Parse GTF
        # Track progress
        message("Parsing GTF...")
        
        # Subset exon entries
        gtf <- gtf[which(gtf$V3=="exon"), ]

        # Retrieve gene names        
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("gene_name", x, value=TRUE))
        . <- gsub("gene_name", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)

        gtf$gene_short_name <- .

        # Convert exon coordinates to intron (to match STAR output)
        gtf$V4 <- gtf$V4 - 1
        gtf$V5 <- gtf$V5 + 1

    # Collapse start position
        # Track progress
        message("Matching gene names with SJ start coordinates in GTF...")
        
        # Keep unique entries
        gtf.small <- gtf[, c("V1", "V5", "gene_short_name")]
        gtf.small$chr.pos <- paste("chr", gtf$V1, ":", gtf$V5, sep="")
        gtf.small <- gtf.small[,c("chr.pos", "gene_short_name")]
        gtf.small <- unique(gtf.small)
        
        # Collapse
        gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
        . <- by(gtf.small[,"gene_short_name"], gtf.small[,"chr.pos"], function(x) {paste(x, collapse="|")})
        gtf.small.collapsed <- data.frame("chr"=as.character(names(.)), "gene_short_name"=as.character(.), stringsAsFactors=FALSE)
        
        # Save as new object
        gtf.small.collapsed.start <- gtf.small.collapsed
        names(gtf.small.collapsed.start) <- paste(names(gtf.small.collapsed.start), ".start", sep="")
        
    # Collapse end position
        # Track progress
        message("Matching gene names with SJ end coordinates in GTF...")
        
        # Keep unique entries
        gtf.small <- gtf[, c("V1", "V4", "gene_short_name")]
        gtf.small$chr.pos <- paste("chr", gtf$V1, ":", gtf$V4, sep="")
        gtf.small <- gtf.small[,c("chr.pos", "gene_short_name")]
        gtf.small <- unique(gtf.small)
        
        # Collapse
        gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
        . <- by(gtf.small[,"gene_short_name"], gtf.small[,"chr.pos"], function(x) {paste(x, collapse="|")})
        gtf.small.collapsed <- data.frame("chr"=as.character(names(.)), "gene_short_name"=as.character(.), stringsAsFactors=FALSE)
        
        # Save as new object
        gtf.small.collapsed.end <- gtf.small.collapsed
        names(gtf.small.collapsed.end) <- paste(names(gtf.small.collapsed.end), ".end", sep="")
        
    # Annotate splice junctions
        # Track progress
        message("Annotating splice junctions...")
        
        # Retrieve SJ start, end coordinates
        df$chr.start <- paste(df$chr, ":", df$start, sep="")
        df$chr.end <- paste(df$chr, ":", df$end, sep="")
        #df$chr.start.end <- paste(df$chr, ":", df$start, ":", df$end, sep="")

        # Indicate id for re-ordering later
        df$id <- c(1:nrow(df))

        # Annotate SJ start,end with gene names
        #dim(gtf.small.collapsed.start) ; dim(gtf.small.collapsed.end) ; dim(df)
        df <- join(df, gtf.small.collapsed.start, by="chr.start", type="left")
        df <- join(df, gtf.small.collapsed.end, by="chr.end", type="left")
        #dim(df) ; sum(is.na(df$gene_short_name.start)) ; sum(is.na(df$gene_short_name.end))

        # Merge gene annotation
            # Create new column
            df$sj.type <- NA
            
            # Single-end record: Same gene at SJ start & end
            index.l <- !grepl("|", df$gene_short_name.start, fixed=TRUE) & !grepl("|", df$gene_short_name.end, fixed=TRUE) & !is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & df$gene_short_name.start==df$gene_short_name.end
            index <- which(index.l==TRUE)
            
            if(length(index) >= 1) {
                
                df$sj.type[index] <- "start_known.single.gene|end_known.single.gene|same"
                
            }
            
            # Single-end record: Different gene at SJ start & end
            index.l <- is.na(df$sj.type) & !grepl("|", df$gene_short_name.start, fixed=TRUE) & !grepl("|", df$gene_short_name.end, fixed=TRUE) & !is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & !df$gene_short_name.start==df$gene_short_name.end
            index <- which(index.l==TRUE)
            
            if(length(index) >= 1) {
                
                df$sj.type[index] <- "start_known.single.gene|end_known.single.gene|different"
                
            }
            
            # No record: Both sj
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & is.na(df$gene_short_name.end)
            index <- which(index.l==TRUE)
            
            if(length(index) >= 1) {
                
                df$sj.type[index] <- "start_unknown.gene|end_unknown.gene"
                
            }
                        
            # No record (start SJ) + end single gene
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & !grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_unknown.gene|end_known.single.gene"
                 
            }
            
            # No record (start SJ) + end multi gene
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_unknown.gene|end_known.multi.gene"
                 
            }
            
            # No record (end SJ) + start single gene
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & !grepl("|", df$gene_short_name.start, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.single.gene|end_unknown.gene"
                 
            }
            
            # No record (end SJ) + start multi gene
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.multi.gene|end_unknown.gene"
                 
            }
            
            # Multiple genes start and end
            index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, fixed=TRUE) & grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.multi.gene|end_known.multi.gene"
                 
            }
            
            # Multiple genes start + single gene end
            index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, fixed=TRUE) & !grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.multi.gene|start_known.single.gene"
                 
            }

            # Multiple genes end + single gene start
            index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & !grepl("|", df$gene_short_name.start, fixed=TRUE) & grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.single.gene|end_known.multi.gene"
                 
            }

            # Check for missing sj annotations
            if(sum(is.na(df$sj.type)) == 0) {
                
                message("All SJ successfully annotated")
                
            } else {
                
                message("Some SJ NOT successfully annotated")
            }

        # Formet input for MARVEL
        df <- df[,c("coord.intron", "gene_short_name.start", "gene_short_name.end", "sj.type")]
    
    #########################################################
    
    # Update slot
    MarvelObject$sj.metadata <- df
            
    # Return final object
    return(MarvelObject)
        
}
