#' @title Plots the locations of specified splice junction relative to isoforms
#'
#' @description Plots the locations of specified splice junction relative to isoforms. List of isoforms are retrieved from GTF.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param coord.intron Character string. Coordinates of splice junction whose splice junction will be plotted.
#' @param coord.intron.ext Numeric value. Number of bases to extend the splice junction start and end coordinates into the exons. Helpful to enhance splice junction locations on the plot. Default is \code{50}.
#' @param rescale_introns Logical value. If set to \code{TRUE}, the intron length will be shorten. Helpful when introns are very long and focus visualisation of exons and splice junctions. Default is \code{FALSE}.
#' @param show.protein.coding.only Logical value. If set to \code{TRUE} (default), only protein-coding isoforms will be displayed.
#' @param anno.label.size Numeric value. Font size of isoform ID labels. Default is \code{3}.
#' @param anno.colors Vector of character strings. Colors for non-coding UTRs, coding exons, and splice junctions, respectively. Default is \code{c("black", "gray", "red")}.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$adhocGene$SJPosition$Plot}, \code{MarvelObject$adhocGene$SJPosition$metadata}, \code{MarvelObject$adhocGene$SJPosition$exonfile}, and \code{MarvelObject$adhocGene$SJPosition$cdsfile}.
#'
#' @importFrom plyr join
#' @importFrom S4Vectors Rle
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom IRanges IRanges
#' @import wiggleplotr
#'
#' @export

adhocGene.PlotSJPosition.10x <- function(MarvelObject, coord.intron, coord.intron.ext=50, rescale_introns=FALSE, show.protein.coding.only=TRUE, anno.label.size=3, anno.colors=c("black", "gray", "red")) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    coord.intron <- coord.intron
    sj.metadata <- MarvelObject$sj.metadata
    gtf <- MarvelObject$gtf
    coord.intron.ext <- coord.intron.ext
    rescale_introns <- rescale_introns
    show.protein.coding.only <- show.protein.coding.only
    anno.label.size <- anno.label.size
    anno.colors <- anno.colors
    
    # Example arguments
    #MarvelObject <- marvel
    #coord.intron <- coord.intron
    #sj.metadata <- MarvelObject$sj.metadata
    #gtf <- MarvelObject$gtf
    #coord.intron.ext <- 10
    #rescale_introns <- TRUE
    #show.protein.coding.only <- FALSE
    #anno.label.size <- 3
    #anno.colors <- c("black", "gray", "red")
    
    #################################################################
    #################### RETRIEVE TRANSCRIPTS #######################
    #################################################################
     
    # Retrieve gene name
    gene_short_name <- sj.metadata[which(sj.metadata$coord.intron==coord.intron), "gene_short_name.start"]
     
    # Subset (by approximation) releavnt gene
    gtf <- gtf[grep(gene_short_name, gtf$V9, fixed=TRUE), ]
    
    # Report progress
    print("Retrieving transcripts from GTF file...")
    
    . <- strsplit(gtf$V9, split=";")
    . <- sapply(., function(x) grep("gene_name", x, value=TRUE))
    . <- gsub("gene_name", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)

    gtf$gene_short_name <- .
    
    # Subset relevant gene
    gtf <- gtf[which(gtf$gene_short_name==gene_short_name), ]
      
    # Retrieve selected attributes
        # gene_id
        #. <- strsplit(gtf$V9, split=";")
        #. <- sapply(., function(x) grep("gene_id", x, value=TRUE))
        #. <- gsub("gene_id", "", .)
        #. <- gsub(" ", "", .)
        #. <- gsub("\"", "", .)

        #gtf$gene_id <- .
        
        # transcript_id
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("transcript_id", x, value=TRUE))
        . <- gsub("transcript_id", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)

        gtf$transcript_id <- .
        
        # transcript_biotype
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("transcript_biotype", x, value=TRUE))
        . <- gsub("transcript_biotype", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)

        if(length(unique(.))==1 & unique(.)[1] == "character(0)") {
                    
            . <- strsplit(gtf$V9, split=";")
            . <- sapply(., function(x) grep("transcript_type", x, value=TRUE))
            . <- gsub("transcript_type", "", .)
            . <- gsub(" ", "", .)
            . <- gsub("\"", "", .)
            
            gtf$transcript_biotype <- .
            
        } else {
            
            gtf$transcript_biotype <- .
            
        }
        
        # exon_id
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("exon_id", x, value=TRUE))
        . <- gsub("exon_id", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)

        gtf$exon_id <- .
    
    #gtf.backup <- gtf
        
    # Annotate transcripts with ORF status
    anno <- unique(gtf[,c("transcript_id", "transcript_biotype")])
    anno <- anno[grep("ENST|ENSMUST", anno$transcript_id), ]
    
    # Report progress
    print(paste(nrow(anno), " transcripts identified", sep=""))
    
    #################################################################
    ############# PREPARE WIGGLEPLOTR INPUT: EXON FILE ##############
    #################################################################
    
    # Define transcript ids
    transcript_ids <- anno$transcript_id
    
    grange.exon.list <- list()
    
    for(i in 1:length(transcript_ids)) {
        
        # Retrieve exons of relevant transcript
        gtf.small <- gtf[which(gtf$transcript_id==transcript_ids[i]), ]
        gtf.small <- gtf.small[which(gtf.small$V3=="exon"), ]
        
        # Create GRange object
        grange <- GRanges(seqnames=Rle(gtf.small$V1),
                          ranges=IRanges(gtf.small$V4,
                                 width=(gtf.small$V5-gtf.small$V4)+1
                                 ),
                          strand=gtf.small$V7[1],
                          exon_id=gtf.small$exon_id,
                          exon_name=gtf.small$exon_id,
                          exon_rank=c(1:length(gtf.small$exon_id))
                          )
                          
        # Save into list
        grange.exon.list[[i]] <- grange
        
    }
        
    # Convert list to GRanges list
    grange.exon.list <- GRangesList(grange.exon.list)
    names(grange.exon.list) <- transcript_ids
    
    #################################################################
    ########### PREPARE WIGGLEPLOTR INPUT: EXON FILE + SJ ###########
    #################################################################
    
    grange.exon.sj.list <- list()
    transcript.ids <- NULL
    
    for(i in 1:length(grange.exon.list)) {
        
        # Retrieve transcript
        grange <- grange.exon.list[[i]]
        exon <- as.data.frame(grange)
        exon$start <- as.numeric(exon$start)
        exon$end <- as.numeric(exon$end)
        
        # Retrieve SJ chr, start, end
        . <- strsplit(coord.intron, split=":", fixed=TRUE)[[1]]
        chr.sj <- .[1]
        start.sj <- as.numeric(.[2])
        end.sj <- as.numeric(.[3])
        
        # Find start, end exon
        exon.sj.start <- exon$end[which(exon$end==start.sj-1)]
        exon.sj.end <- exon$start[which(exon$start==end.sj+1)]
        
        # Retrieve IRanges to subset
        exon.small <- exon[which(exon$end==exon.sj.start | exon$start==exon.sj.end),  ]
        
        # Trim exon length
        exon.small$start[which(exon.small$end==exon.sj.start)] <- exon.small$end[which(exon.small$end==exon.sj.start)] - coord.intron.ext
        exon.small$end[which(exon.small$start==exon.sj.end)] <- exon.small$start[which(exon.small$start==exon.sj.end)] + coord.intron.ext
        
        if(nrow(exon.small) != 0) {
            
            # Create GRanges object
            grange.sj <- GRanges(seqnames=Rle(exon.small$seqnames),
                              ranges=IRanges(exon.small$start,
                                     width=(exon.small$end-exon.small$start)+1
                                     ),
                              strand=exon.small$strand,
                              exon_id=exon.small$exon_id,
                              exon_name=exon.small$exon_id,
                              exon_rank=c(1:length(exon.small$exon_id))
                              )
            
            
            grange.exon.sj.list[[i]] <- grange.sj
            transcript.ids[i] <- names(grange.exon.list)[i]
                                          
        } else {
            
            grange.exon.sj.list[[i]] <- FALSE
            transcript.ids[i] <- FALSE
            
        }
        
    }
        
    # Convert list to GRanges list
    index.keep <- which(transcript.ids != FALSE)
    grange.exon.sj.list <- grange.exon.sj.list[index.keep]
    transcript.ids <- transcript.ids[index.keep]
    
    if(length(grange.exon.sj.list) != 0) {
        
        grange.exon.sj.list <- GRangesList(grange.exon.sj.list)
        names(grange.exon.sj.list) <- transcript.ids
    
    }
    
    #################################################################
    ################ PREPARE WIGGLEPLOTR INPUT: CDS #################
    #################################################################
    
    # Define transcript ids
    transcript_ids <- anno$transcript_id
    
    grange.cds.list <- list()
    transcript.ids <- NULL
    
    for(i in 1:length(transcript_ids)) {
        
        # Retrieve exons of relevant transcript
        gtf.small <- gtf[which(gtf$transcript_id==transcript_ids[i]), ]
        gtf.small <- gtf.small[which(gtf.small$V3=="CDS"), ]
        
        if(nrow(gtf.small) != 0) {
            
            # Create GRange object
            grange <- GRanges(seqnames=Rle(gtf.small$V1),
                              ranges=IRanges(gtf.small$V4,
                                     width=(gtf.small$V5-gtf.small$V4)+1
                                     ),
                              strand=gtf.small$V7[1],
                              exon_id=gtf.small$exon_id,
                              exon_name=gtf.small$exon_id,
                              exon_rank=c(1:length(gtf.small$exon_id))
                              )
                              
            # Save into list
            grange.cds.list[[i]] <- grange
            transcript.ids[i] <- transcript_ids[i]
        
        } else {
            
            # Save into list
            grange.cds.list[[i]] <- FALSE
            transcript.ids[i] <- FALSE
            
        }
        
    }
        
    # Convert list to GRanges list
    index.keep <- which(transcript.ids != FALSE)
    grange.cds.list <- grange.cds.list[index.keep]
    transcript.ids <- transcript.ids[index.keep]
    
    if(length(grange.cds.list) != 0) {
        
        grange.cds.list <- GRangesList(grange.cds.list)
        names(grange.cds.list) <- transcript.ids
    
    }
    
    #################################################################
    ############# PREPARE WIGGLEPLOTR INPUT: METADATA ###############
    #################################################################
    
    # Retrieve transcript ids
    metadata <- data.frame("transcript_id"=names(grange.exon.list), stringsAsFactors=FALSE)

    # Annotate gene id, name
    #metadata$gene_id <- gene_id
    metadata$gene_short_name <- gene_short_name
    
    # Indicate strand
    strand <- gtf$V7[1]
    
    if(strand=="+") {
        
        metadata$strand <- 1
        
    } else {
        
        metadata$strand <- -1
    
    }
    
    # Annotate transcript type
        # Retrieve transcript type
        metadata <- plyr::join(metadata, anno, by="transcript_id", type="left")
        metadata$transcript_id.biotype <- paste(metadata$transcript_id, " (", metadata$transcript_biotype, ")", sep="")
        
        # Update GRange exon list
        . <- data.frame("transcript_id"=names(grange.exon.list), stringsAsFactors=FALSE)
        . <- plyr::join(., metadata[,c("transcript_id", "transcript_id.biotype")], by="transcript_id", type="left")
        names(grange.exon.list) <- .$transcript_id.biotype
        
        # Update GRange exon-SJ list
        if(length(grange.exon.sj.list) != 0) {
            
            . <- data.frame("transcript_id"=names(grange.exon.sj.list), stringsAsFactors=FALSE)
            . <- plyr::join(., metadata[,c("transcript_id", "transcript_id.biotype")], by="transcript_id", type="left")
            names(grange.exon.sj.list) <- paste(.$transcript_id.biotype, "_SJ", sep="")
        
        }
        
        # Update GRange CDS list
        . <- data.frame("transcript_id"=names(grange.cds.list), stringsAsFactors=FALSE)
        . <- plyr::join(., metadata[,c("transcript_id", "transcript_id.biotype")], by="transcript_id", type="left")
        names(grange.cds.list) <- .$transcript_id.biotype
        
        # Update metadata
        metadata$transcript_id <- metadata$transcript_id.biotype
        
        # Remove intermediate column
        metadata$transcript_id.biotype <- NULL
        metadata$transcript_biotype <- NULL
        
    # Update exon file, metadata with exon-SJ
        # Update exon file
        if(length(grange.exon.sj.list) != 0) {
            
            grange.exon.list <- c(grange.exon.list, grange.exon.sj.list)
            
        }
        
        # Update metadata file
        if(length(grange.exon.sj.list) != 0) {
            
            #metadata <- data.frame("transcript_id"=names(grange.exon.list),
                            #"gene_id"=gene_id,
                            #"gene_short_name"=gene_short_name,
                            #stringsAsFactors=FALSE
                            #)
                            
            metadata <- data.frame("transcript_id"=names(grange.exon.list),
                            "gene_short_name"=gene_short_name,
                            stringsAsFactors=FALSE
                            )
                            
            strand <- gtf$V7[1]
            
            if(strand=="+") {
                
                metadata$strand <- 1
                
            } else {
                
                metadata$strand <- -1
            
            }
            
        }
    
    #################################################################
    ################# FILTER FOR SPECIFIC BIOTYPE ###################
    #################################################################
    
    if(show.protein.coding.only==TRUE) {
        
        transcript_ids <- metadata[grep("protein_coding", metadata$transcript_id, fixed=TRUE), "transcript_id"]
        
        if(length(transcript_ids) != 0) {
            
            metadata <- metadata[grep("protein_coding", metadata$transcript_id, fixed=TRUE), ]
            grange.exon.list <- grange.exon.list[metadata$transcript_id]
            
            overlap <- intersect(names(grange.cds.list), transcript_ids)
            grange.cds.list <- grange.cds.list[overlap]
            
            
            
        } else {
            
            print("No protein-coding transcripts found for this gene")
            MarvelObject$adhocGene$SJPosition$metadata <- metadata
            return
            
        }
        
        
    }
    
    #################################################################
    ######################## WIGGLEPLOTR ############################
    #################################################################
  
    plot <- plotTranscripts(exons=grange.exon.list,
                            cdss=grange.cds.list,
                            transcript_annotations=metadata,
                            rescale_introns=rescale_introns,
                            new_intron_length=50,
                            flanking_length=c(50,50),
                            connect_exons=TRUE,
                            transcript_label=TRUE,
                            region_coords = NULL,
                            anno.colors=anno.colors,
                            anno.label.size=anno.label.size
                            )
                            
    # Save into new slots
    MarvelObject$adhocGene$SJPosition$Plot <- plot
    MarvelObject$adhocGene$SJPosition$metadata <- metadata
    MarvelObject$adhocGene$SJPosition$exonfile <- grange.exon.list
    MarvelObject$adhocGene$SJPosition$cdsfile <- grange.cds.list
    
    # Return final object
    return(MarvelObject)
            
}


