#' @title Find premature terminal codon (PTC) for retained-intron (RI) located on the negative strand of the transcript
#'
#' @description Finds PTC(s) introduced by alternative exons into protein-coding transcripts.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues.PSI} and \code{ParseGTF} function.
#' @param tran_id Character string. Vector of \code{tran_id} to look for PTCs.
#' @param gene_id Character string. Vector of \code{gene_id} corresponding to the \code{tran_id} argument.
#'
#' @return A data frame of transcripts containing splicing events meeting the \code{psi.de.sig} and \code{psi.de.diff} criteria are categorised based on the presence or absence of PTCs.
#'
#' @importFrom plyr join
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' # Define relevant event type
#' results <- marvel.demo$DE$PSI$Table[["ad"]]
#' index.1 <- which(results$event_type=="RI")
#' index.2 <- grep(":-@", results$tran_id, fixed=TRUE)
#' index <- intersect(index.1, index.2)
#' results <- results[index, ]
#' tran_id <- results$tran_id[1]
#' gene_id <- results$gene_id[1]
#'
#' # Find PTC
#' results <- FindPTC.RI.NegStrand(MarvelObject=marvel.demo,
#'                                 tran_id=tran_id,
#'                                 gene_id=gene_id
#'                                 )
#'
#' # Check output
#' head(results)

FindPTC.RI.NegStrand <- function(MarvelObject, tran_id, gene_id) {

    # Define arguments
    gtf <- MarvelObject$NMD$GTF
    tran_id <- tran_id
    gene_id <- gene_id

    # Example arguments
    #gtf <- gtf
    #tran_id <- tran_ids[3]
    #gene_id <- gene_ids[3]

    # Create container to keep results
    results.list <- list()
    
    ##################################################################################
    ################################# DEFINE EXONS ###################################
    ##################################################################################

    # Retrieve chr
    . <- strsplit(tran_id, split=":-@", fixed=TRUE)[[1]][2]
    chr <- strsplit(., split=":", fixed=TRUE)[[1]][1]

    # 5' cons'
    . <- strsplit(tran_id, split=":-@", fixed=TRUE)[[1]][1]
    cons.exon.5.start <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][2])
    cons.exon.5.end <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][3])

    # 3' cons'
    . <- strsplit(tran_id, split=":-@", fixed=TRUE)[[1]][2]
    cons.exon.3.start <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][2])
    cons.exon.3.end <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][3])

    # alt exon
    . <- strsplit(tran_id, split=":-@", fixed=TRUE)[[1]][1]
    alt.exon.start <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][3]) - 1

    . <- strsplit(tran_id, split=":-@", fixed=TRUE)[[1]][2]
    alt.exon.end <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][2]) + 1

    ##################################################################################
    ##################### RETRIEVE TRANSCRIPTS WITHOUT ALT. EXON #####################
    ##################################################################################

    # Subset relevant gene
    gtf.small <- gtf[which(gtf$gene_id==gene_id), ]

    # Subset relevant transcripts
        # Retrieve 5' cons. exon
        index.5 <- which(gtf.small$V4==cons.exon.5.end)

        # Check and retrieve adjacent exon
        adjacent.start <- gtf.small$V5[index.5 + 1]
        index.adjacent.start <- which(adjacent.start==cons.exon.3.start)

        # Subset relevant transcripts + track progress
        index.5 <- index.5[index.adjacent.start]
        transcript_ids <- gtf.small[index.5, "transcript_id"]
        gtf.small <- gtf.small[which(gtf.small$transcript_id %in% transcript_ids), ]

        if(nrow(gtf.small) == 0) {
            
            #message("No transcripts with matching SJ found for this event")
            
            results <- data.frame("tran_id"=tran_id,
                                  "gene_id"=gene_id,
                                  "transcript_id"=NA,
                                  "aa.length"=NA,
                                  "stop.position.aa"=NA,
                                  "sj.last.position.aa"=NA,
                                  "ptc.sj.last.distance"=NA,
                                  "NMD"="No transcripts with matching SJ",
                                  "event_type"="RI",
                                  "strand"="-",
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
                                            "event_type"="RI",
                                            "strand"="-",
                                            stringsAsFactors=FALSE
                                            )
            
        }
        
        # Subset relevant transcripts
        gtf.small <- gtf.small[which(gtf.small$transcript_type=="protein_coding"), ]

        if(nrow(gtf.small) == 0) {
            
            #message("No protein-coding transcripts found for this gene")
            
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
                                            "event_type"="RI",
                                            "strand"="-",
                                            stringsAsFactors=FALSE
                                            )
            
        }

        # Subset relevant transcripts
        transcript_ids <- transcript_ids[index.keep==TRUE]

        gtf.small <- gtf.small[which(gtf.small$transcript_id %in% transcript_ids), ]

        if(nrow(gtf.small) == 0) {
            
            #message("No transcripts with both start and stop codon found for this gene")
            
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
        index.5 <- cons.exon.5.end <= gtf.small.$V4[1]

        # Check if 3' cons. exon comes before CDS
        index.3 <- cons.exon.3.start >= gtf.small.$V5[nrow(gtf.small.)]
        
        # Indicate in or out of ORF
        if(index.5==TRUE & index.3==TRUE) {
        
            index.keep[i] <- TRUE
        
        } else {
        
            index.keep[i] <- FALSE
        
        }

    }

    # Record irrelevant transcripts
    transcript_ids.null <- transcript_ids[-index.keep==FALSE]

    if(length(transcript_ids.null) != 0) {
        
        results.list[[3]] <- data.frame("tran_id"=tran_id,
                                        "gene_id"=gene_id,
                                        "transcript_id"=transcript_ids.null,
                                        "aa.length"=NA,
                                        "stop.position.aa"=NA,
                                        "sj.last.position.aa"=NA,
                                        "ptc.sj.last.distance"=NA,
                                        "NMD"="splicing event located outside of ORF",
                                        "event_type"="RI",
                                        "strand"="-",
                                        stringsAsFactors=FALSE
                                        )
        
    }
    
    # Subset relevant transcripts
    transcript_ids <- transcript_ids[index.keep==TRUE]

    gtf.small <- gtf.small[which(gtf.small$transcript_id %in% transcript_ids), ]
    #message(paste(length(transcript_ids), " relevant transcripts identified", sep=""))

    if(nrow(gtf.small) == 0) {
        
        #message("No transcripts with both start and stop codon found for this gene")
        
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
        index <- which(gtf.small.$V4==cons.exon.5.end)
        gtf.small.5 <- gtf.small.[c(1:index), ]
        
        # Retrieve 3' cons. exon up to and including last exon
        index <- which(gtf.small.$V5==cons.exon.3.start)
        gtf.small.3 <- gtf.small.[c(index:nrow(gtf.small.)), ]
        
        # Create new row for alt exon
        gtf.small.alt <- data.frame("V1"=chr,
                                    "V4"=alt.exon.end,
                                    "V5"=alt.exon.start,
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
        DNAString <- BSgenome::getSeq(BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens, chr, start=start, end=end)
        DNAString <- Biostrings::reverseComplement(DNAString) # Only for -ve strand!
        DNAString <- paste(as.character(DNAString), collapse="")
        
        # Translate sequence
        dna <- Biostrings::DNAStringSet(DNAString)
        aa <- Biostrings::translate(dna)
        amino.acids[i] <- as.character(aa)
        
        # Estimate position of STOP codon
        stop.position.aa[i] <- gregexpr(pattern="*", aa, fixed=TRUE)[[1]][1]
                                
        # Retrieve position of last SJ
        sj.last.position.aa[i] <- nchar(amino.acids[i]) -
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
    results$event_type <- "RI"
    results$strand <- "-"
    
    # Save results
    results.list[[4]] <- results

    # Merge all result tables
    results.final <- do.call(rbind.data.frame, results.list)
    
    # Return MARVEL object
    return(results.final)

}
