#' @title Find Premature Terminal Codon (PTC) for Alternative 5' Spice Site (A5SS) Located on the Positive Strand of the Transcript
#'
#' @description
#' \code{FindPTC.A5SS.PosStrand} finds PTC(s) introduced by alternative exons into protein-coding transcripts.
#'
#' @details
#' This function finds PTC(s) introduced by alternative exons into protein-coding transcripts. It also records the distance between a PTCs and the final splice junction for a given protein-coding transcript. Non-protein-coding transcripts or transcripts in which splicing events are located outside of the transcripts' open-reading frame (ORF) are not analysed for PTCs but are noted.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues.PSI} and \code{ParseGTF} function.
#' @param tran_id Character string. Vector of \code{tran_id} to look for PTCs.
#' @param gene_id Character string. Vector of \code{gene_id} corresponding to the \code{tran_id} argument.
#' @export
#'
#' @return A data frame of transcripts containing splicing events meeting the \code{psi.de.sig} and \code{psi.de.diff} criteria are categorised based on the presence or absence of PTCs.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @import BSgenome.Hsapiens.NCBI.GRCh38
#' @importFrom Biostrings DNAStringSet translate reverseComplement
#' @importFrom BSgenome getSeq
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Retrieve features
#' df <- marvel$DE$PSI$Table
#' df <- df[which(df$event_type=="A5SS"), ]
#' df <- df[grep(":+@", df$tran_id, fixed=TRUE), ]
#' tran_id <- df$tran_id[1]
#' gene_id <- df$gene_id[1]
#'
#' # Run example
#' results <- FindPTC.A5SS.PosStrand(MarvelObject=marvel,
#'                                  tran_id=tran_id,
#'                                  gene_id=gene_id
#'                                  )
#'
#' # Check output
#' results

FindPTC.A5SS.PosStrand <- function(MarvelObject, tran_id, gene_id) {

    # Define arguments
    gtf <- MarvelObject$NMD$GTF
    tran_id <- tran_id
    gene_id <- gene_id

    # Example arguments
    #gtf <- gtf
    #tran_id <- tran_ids[1]
    #gene_id <- gene_ids[1]
    
    # Create container to keep results
    results.list <- list()
    
    ##################################################################################
    ################################# DEFINE EXONS ###################################
    ##################################################################################

    # Retrieve chr
    . <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][2]
    chr <- strsplit(., split=":", fixed=TRUE)[[1]][1]

    # 5' cons'
    . <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][1]
    . <- strsplit(., split="|", fixed=TRUE)[[1]][1]
    cons.exon.5.start <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][2])
    cons.exon.5.end <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][3])

    # 3' cons'
    . <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][2]
    cons.exon.3.start <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][2])
    cons.exon.3.end <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][3])

    # alt exon
    alt.exon.start <- cons.exon.5.end + 1

    . <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][1]
    alt.exon.end <- as.numeric(strsplit(., split="|", fixed=TRUE)[[1]][2])

    ##################################################################################
    ##################### RETRIEVE TRANSCRIPTS WITHOUT ALT. EXON #####################
    ##################################################################################

    # Subset relevant gene
    gtf.small <- gtf[which(gtf$gene_id==gene_id), ]

    # Subset relevant transcripts
        # Retrieve 5' cons. exon
        index.5 <- which(gtf.small$V5==cons.exon.5.end)

        # Check and retrieve adjacent exon
        adjacent.start <- gtf.small$V4[index.5 + 1]
        index.adjacent.start <- which(adjacent.start==cons.exon.3.start)

        # Subset relevant transcripts + track progress
        index.5 <- index.5[index.adjacent.start]
        transcript_ids <- gtf.small[index.5, "transcript_id"]
        gtf.small <- gtf.small[which(gtf.small$transcript_id %in% transcript_ids), ]

        if(nrow(gtf.small) == 0) {
            
            #print("No transcripts with matching SJ found for this event")
            
            results <- data.frame("tran_id"=tran_id,
                                  "gene_id"=gene_id,
                                  "transcript_id"=NA,
                                  "aa.length"=NA,
                                  "stop.position.aa"=NA,
                                  "sj.last.position.aa"=NA,
                                  "ptc.sj.last.distance"=NA,
                                  "NMD"="No transcripts with matching SJ",
                                  "event_type"="A5SS",
                                  "strand"="+",
                                  stringsAsFactors=FALSE
                                  )
            
            return(results)
        
        }
        
    # Subset protein-coding transcripts
        # Record irrelevant transcripts
        gtf.small.null <- gtf.small[which(gtf.small$transcript_type != "protein_coding"), ]

        if(nrow(gtf.small.null) != 0) {
            
            results.list[[1]] <- data.frame("tran_id"=tran_id,
                                            "gene_id"=gene_id,
                                            "transcript_id"=unique(gtf.small.null$transcript_id),
                                            "aa.length"=NA,
                                            "stop.position.aa"=NA,
                                            "sj.last.position.aa"=NA,
                                            "ptc.sj.last.distance"=NA,
                                            "NMD"="non-protein-coding transcript",
                                            "event_type"="A5SS",
                                            "strand"="+",
                                            stringsAsFactors=FALSE
                                            )
            
        }
        
        # Subset relevant transcripts
        gtf.small <- gtf.small[which(gtf.small$transcript_type=="protein_coding"), ]
        
        if(nrow(gtf.small) == 0) {
            
            #print("No protein-coding transcripts found for this gene")
            
            return(do.call(rbind.data.frame, results.list))
        
        }

    # Subset transcripts with known start and stop codon
        # Check transcripts
        transcript_ids <- unique(gtf.small$transcript_id)

        index.keep <- NULL

        for(i in 1:length(transcript_ids)) {

            # Subset relevant transcript
            gtf.small. <- gtf.small[which(gtf.small$transcript_id==transcript_ids[i]), ]
            
            # Check for start codon
            codon.start <- which(gtf.small.$V3=="start_codon")
            
            # Check for stop codon
            codon.stop <- which(gtf.small.$V3=="stop_codon")
            
            # Indicate if both codons found
            if(length(codon.start)==1 & length(codon.stop)==1) {
            
            
                index.keep[i] <- TRUE
            
            } else {
            
                index.keep[i] <- FALSE
                    
            }

        }

        # Record irrelevant transcripts
        transcript_ids.null <- transcript_ids[index.keep==FALSE]
        
        if(length(transcript_ids.null) != 0) {
            
            results.list[[2]] <- data.frame("tran_id"=tran_id,
                                            "gene_id"=gene_id,
                                            "transcript_id"=transcript_ids.null,
                                            "aa.length"=NA,
                                            "stop.position.aa"=NA,
                                            "sj.last.position.aa"=NA,
                                            "ptc.sj.last.distance"=NA,
                                            "NMD"="transcript with no START and/or STOP codon",
                                            "event_type"="A5SS",
                                            "strand"="+",
                                            stringsAsFactors=FALSE
                                            )
            
        }
        
        # Subset relevant transcripts
        transcript_ids <- transcript_ids[index.keep==TRUE]

        gtf.small <- gtf.small[which(gtf.small$transcript_id %in% transcript_ids), ]
        
        if(nrow(gtf.small) == 0) {
            
            #print("No transcripts with both start and stop codon found for this gene")
            
            return(do.call(rbind.data.frame, results.list))
        
        }

    # Subset CDS entries only
    gtf.small <- gtf.small[which(gtf.small$V3=="CDS"), ]

    # Keep relevant columns
    gtf.small <- gtf.small[,c("V1", "V4", "V5", "gene_id", "transcript_id", "transcript_type")]

    ##################################################################################
    ############### SUBSET TRANSCRIPTS: ALT EXON LOCATED INSIDE ORF ##################
    ##################################################################################

    index.keep <- NULL

    for(i in 1:length(transcript_ids)) {

        # Subset relevant transcript
        gtf.small. <- gtf.small[which(gtf.small$transcript_id==transcript_ids[i]), ]
        
        # Check if 5' cons. exon comes after CDS
        index.5 <- cons.exon.5.end >= gtf.small.$V5[1]

        # Check if 3' cons. exon comes before CDS
        index.3 <- cons.exon.3.start <= gtf.small.$V4[nrow(gtf.small.)]
        
        # Indicate in or out of ORF
        if(index.5==TRUE & index.3==TRUE) {
        
            index.keep[i] <- TRUE
        
        } else {
            
            index.keep[i] <- FALSE
            
        }

    }

    # Record irrelevant transcripts
    transcript_ids.null <- transcript_ids[index.keep==FALSE]

    if(length(transcript_ids.null) != 0) {
        
        results.list[[3]] <- data.frame("tran_id"=tran_id,
                                        "gene_id"=gene_id,
                                        "transcript_id"=transcript_ids.null,
                                        "aa.length"=NA,
                                        "stop.position.aa"=NA,
                                        "sj.last.position.aa"=NA,
                                        "ptc.sj.last.distance"=NA,
                                        "NMD"="splicing event located outside of ORF",
                                        "event_type"="A5SS",
                                        "strand"="+",
                                        stringsAsFactors=FALSE
                                        )
                                        
        
    }

    # Subset relevant transcripts
    transcript_ids <- transcript_ids[index.keep==TRUE]

    gtf.small <- gtf.small[which(gtf.small$transcript_id %in% transcript_ids), ]
    
    #print(paste(length(transcript_ids), " relevant transcripts identified", sep=""))

    if(nrow(gtf.small) == 0) {
        
        #print("No transcripts with both start and stop codon found for this gene")
        
        return(do.call(rbind.data.frame, results.list))
    
    }

    ##################################################################################
    ############################## INSERT ALT EXON ###################################
    ##################################################################################

    .list <- list()

    for(i in 1:length(transcript_ids)) {

        # Subset relevant transcript
        gtf.small. <- gtf.small[which(gtf.small$transcript_id==transcript_ids[i]), ]
        
        # Retrieve exons up to and including 5' cons. exon
        index <- which(gtf.small.$V5==cons.exon.5.end)
        gtf.small.5 <- gtf.small.[c(1:index), ]
        
        # Retrieve 3' cons. exon up to and including last exon
        index <- which(gtf.small.$V4==cons.exon.3.start)
        gtf.small.3 <- gtf.small.[c(index:nrow(gtf.small.)), ]
        
        # Create new row for alt exon
        gtf.small.alt <- data.frame("V1"=chr,
                                    "V4"=alt.exon.start,
                                    "V5"=alt.exon.end,
                                    "gene_id"=gtf.small.5$gene_id[1],
                                    "transcript_id"=gtf.small.5$transcript_id[1],
                                    "transcript_type"=gtf.small.5$transcript_type[1],
                                    stringsAsFactors=FALSE
                                    )
        # Merge
        .list[[i]] <- rbind.data.frame(gtf.small.5, gtf.small.alt, gtf.small.3)

    }

    gtf.small <- do.call(rbind.data.frame, .list)

    ##################################################################################
    ############################# TRANSLATE CDS ######################################
    ##################################################################################

    # Translate CDS and find PTC
    amino.acids <- NULL
    stop.position.aa <- NULL
    sj.last.position.aa <- NULL

    for(i in 1:length(transcript_ids)) {

        # Subset relevant transcript
        gtf.small. <- gtf.small[which(gtf.small$transcript_id==transcript_ids[i]), ]
        
        # Define chr, start, end
        chr <- gsub("chr", "", gtf.small.$V1, fixed=TRUE)
        start <-  gtf.small.$V4
        end <-  gtf.small.$V5
        
        # Retrieve nt sequence
        DNAString <- getSeq(Hsapiens, chr, start=start, end=end)
        DNAString <- paste(as.character(DNAString), collapse="")
        
        # Translate sequence
        dna <- DNAStringSet(DNAString)
        aa <- translate(dna)
        amino.acids[i] <- as.character(aa)
        
        # Retrieve position of 1st STOP codon
        stop.position.aa[i] <- gregexpr(pattern="*", aa, fixed=TRUE)[[1]][1]
        
        # Retrieve position of last SJ
        sj.last.position.aa[i] <-   nchar(amino.acids[i]) -
                                    floor( (gtf.small.$V5[nrow(gtf.small.)] - gtf.small.$V4[nrow(gtf.small.)]) / 3)
     
    }

    # Tabulate
    results <- data.frame("tran_id"=tran_id,
                          "gene_id"=gene_id,
                          "transcript_id"=transcript_ids,
                          "aa.length"=nchar(amino.acids),
                          "stop.position.aa"=stop.position.aa,
                          "sj.last.position.aa"=sj.last.position.aa,
                          stringsAsFactors=FALSE
                          )
                          
    # Compute relative position of PTC to last SJ
    results$ptc.sj.last.distance <- (results$sj.last.position.aa - results$stop.position.aa) * 3

    # Predict NMD
    results$NMD <- ifelse(results$ptc.sj.last.distance > 50,
                          "yes_distance of PTC to final SJ > 50bp",
                          "no_distance of PTC to final SJ <= 50bp"
                          )

    # Recode selected variables when no PTC identified
    results$ptc.sj.last.distance[which(results$stop.position.aa==-1)] <- NA
    results$NMD[which(results$stop.position.aa==-1)] <- "no_PTC not found"

    # Annotate event type and strand
    results$event_type <- "A5SS"
    results$strand <- "+"
    
    # Save results
    results.list[[4]] <- results

    # Merge all result tables
    results.final <- do.call(rbind.data.frame, results.list)

    # Return MARVEL object
    return(results.final)
    
    # Return MARVEL object
    return(results)

}
