#' @title Find Premature Terminal Codon (PTC)
#'
#' @description
#' \code{FindPTC} finds PTC(s) introduced by alternative exons into protein-coding transcripts.
#'
#' @details
#' This function finds PTC(s) introduced by alternative exons into protein-coding transcripts. It also records the distance between a PTCs and the final splice junction for a given protein-coding transcript. Non-protein-coding transcripts or transcripts in which splicing events are located outside of the transcripts' open-reading frame (ORF) are not analysed for PTCs but are noted.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues.PSI} and \code{ParseGTF} functions.
#' @param psi.de.sig Numeric value. Adjusted p-value below which the splicing event will be analysed for PTCs.
#' @param psi.de.diff Numeric value. Positive delta percent spliced-in (PSI) value above which the splicing event will be analysed for PTCs. "Positive" because only an increased in PSI value leads to increased alternative exon inclusion in the transcript.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$NMD$Prediction}. Transcripts containing splicing events meeting the \code{psi.de.sig} and \code{psi.de.diff} criteria are categorised based on the presence or absence of PTCs.
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
#' # Run example
#' marvel <- FindPTC(MarvelObject=marvel,
#'                   psi.de.sig=0.10,
#'                   psi.de.diff=0.05
#'                   )
#'
#' # Check output
#' marvel$NMD$Prediction


FindPTC <- function(MarvelObject, psi.de.sig, psi.de.diff) {

    # Define arguments
    de.psi <- MarvelObject$DE$PSI$Table
    psi.de.sig <- psi.de.sig
    psi.de.diff <- psi.de.diff

    # Example arguments
    #MarvelObject <- marvel
    #de.psi <- marvel$DE$PSI$Table
    #psi.de.sig <- 0.10
    #psi.de.diff <- 0.05
    
    # Subset sig. events
    de.psi.small <- de.psi[which(de.psi$p.val.adj < psi.de.sig & de.psi$mean.diff > psi.de.diff), ]
    
    # Remove MXE events
    de.psi.small <- de.psi.small[which(de.psi.small$event_type != "MXE"), ]
    
    # Create container to keep results
    results.list <- list()

    #####################################################################

    # Predict NMD: SE +ve strand
    index.1 <- which(de.psi.small$event_type=="SE")
    index.2 <- grep(":+@", de.psi.small$tran_id, fixed=TRUE)
    index <- intersect(index.1, index.2)

    tran_ids <- de.psi.small$tran_id[index]
    gene_ids <- de.psi.small$gene_id[index]
    
    print(paste(length(tran_ids),  " SE +ve strand identified", sep=""))
    
    if(length(tran_ids) != 0) {
        
        .list <- list()
            
        for(i in 1:length(tran_ids)) {

            .list[[i]] <- FindPTC.SE.PosStrand(MarvelObject, tran_id=tran_ids[i], gene_ids[i])
            
        }
        
        results.list[[1]] <- do.call(rbind.data.frame, .list)
    
    }
    
    # Predict NMD: SE -ve strand
    index.1 <- which(de.psi.small$event_type=="SE")
    index.2 <- grep(":-@", de.psi.small$tran_id, fixed=TRUE)
    index <- intersect(index.1, index.2)

    tran_ids <- de.psi.small$tran_id[index]
    gene_ids <- de.psi.small$gene_id[index]
    
    print(paste(length(tran_ids),  " SE -ve strand identified", sep=""))
    
    if(length(tran_ids) != 0) {
        
        .list <- list()
            
        for(i in 1:length(tran_ids)) {

            .list[[i]] <- FindPTC.SE.NegStrand(MarvelObject, tran_id=tran_ids[i], gene_ids[i])
            
        }
        
        results.list[[2]] <- do.call(rbind.data.frame, .list)
    
    }
    
    # Predict NMD: RI +ve strand
    index.1 <- which(de.psi.small$event_type=="RI")
    index.2 <- grep(":+@", de.psi.small$tran_id, fixed=TRUE)
    index <- intersect(index.1, index.2)

    tran_ids <- de.psi.small$tran_id[index]
    gene_ids <- de.psi.small$gene_id[index]
    
    print(paste(length(tran_ids),  " RI +ve strand identified", sep=""))
    
    if(length(tran_ids) != 0) {
        
        .list <- list()
            
        for(i in 1:length(tran_ids)) {

            .list[[i]] <- FindPTC.RI.PosStrand(MarvelObject, tran_id=tran_ids[i], gene_ids[i])
            
        }
        
        results.list[[3]] <- do.call(rbind.data.frame, .list)
    
    }
    
    # Predict NMD: RI -ve strand
    index.1 <- which(de.psi.small$event_type=="RI")
    index.2 <- grep(":-@", de.psi.small$tran_id, fixed=TRUE)
    index <- intersect(index.1, index.2)

    tran_ids <- de.psi.small$tran_id[index]
    gene_ids <- de.psi.small$gene_id[index]
    
    print(paste(length(tran_ids),  " RI -ve strand identified", sep=""))
    
    if(length(tran_ids) != 0) {
        
        .list <- list()
            
        for(i in 1:length(tran_ids)) {

            .list[[i]] <- FindPTC.RI.NegStrand(MarvelObject, tran_id=tran_ids[i], gene_ids[i])
            
        }
        
        results.list[[4]] <- do.call(rbind.data.frame, .list)
    
    }
    
    # Predict NMD: A5SS +ve strand
    index.1 <- which(de.psi.small$event_type=="A5SS")
    index.2 <- grep(":+@", de.psi.small$tran_id, fixed=TRUE)
    index <- intersect(index.1, index.2)

    tran_ids <- de.psi.small$tran_id[index]
    gene_ids <- de.psi.small$gene_id[index]
    
    print(paste(length(tran_ids),  " A5SS +ve strand identified", sep=""))

    if(length(tran_ids) != 0) {
        
        .list <- list()
            
        for(i in 1:length(tran_ids)) {

            .list[[i]] <- FindPTC.A5SS.PosStrand(MarvelObject, tran_id=tran_ids[i], gene_ids[i])
            
        }
        
        results.list[[5]] <- do.call(rbind.data.frame, .list)
    
    }
    
    # Predict NMD: A5SS -ve strand
    index.1 <- which(de.psi.small$event_type=="A5SS")
    index.2 <- grep(":-@", de.psi.small$tran_id, fixed=TRUE)
    index <- intersect(index.1, index.2)

    tran_ids <- de.psi.small$tran_id[index]
    gene_ids <- de.psi.small$gene_id[index]
    
    print(paste(length(tran_ids),  " A5SS -ve strand identified", sep=""))

    if(length(tran_ids) != 0) {
        
        .list <- list()
            
        for(i in 1:length(tran_ids)) {

            .list[[i]] <- FindPTC.A5SS.NegStrand(MarvelObject, tran_id=tran_ids[i], gene_ids[i])
            
        }
        
        results.list[[6]] <- do.call(rbind.data.frame, .list)
    
    }

    # Predict NMD: A3SS +ve strand
    index.1 <- which(de.psi.small$event_type=="A3SS")
    index.2 <- grep(":+@", de.psi.small$tran_id, fixed=TRUE)
    index <- intersect(index.1, index.2)

    tran_ids <- de.psi.small$tran_id[index]
    gene_ids <- de.psi.small$gene_id[index]
    
    print(paste(length(tran_ids),  " A3SS +ve strand identified", sep=""))

    if(length(tran_ids) != 0) {
        
        .list <- list()
            
        for(i in 1:length(tran_ids)) {

            .list[[i]] <- FindPTC.A3SS.PosStrand(MarvelObject, tran_id=tran_ids[i], gene_ids[i])
            
        }
        
        results.list[[7]] <- do.call(rbind.data.frame, .list)
    
    }
    
    # Predict NMD: A3SS -ve strand
    index.1 <- which(de.psi.small$event_type=="A3SS")
    index.2 <- grep(":-@", de.psi.small$tran_id, fixed=TRUE)
    index <- intersect(index.1, index.2)

    tran_ids <- de.psi.small$tran_id[index]
    gene_ids <- de.psi.small$gene_id[index]
    
    print(paste(length(tran_ids),  " A3SS -ve strand identified", sep=""))

    if(length(tran_ids) != 0) {
        
        .list <- list()
            
        for(i in 1:length(tran_ids)) {

            .list[[i]] <- FindPTC.A3SS.NegStrand(MarvelObject, tran_id=tran_ids[i], gene_ids[i])
            
        }
        
        results.list[[8]] <- do.call(rbind.data.frame, .list)
    
    }

    ######################################################################

    # Merge results
    results <- do.call(rbind.data.frame, results.list)

    # Save to new slot
    MarvelObject$NMD$Prediction <- results
    
    # Return MARVEL object
    return(MarvelObject)
    
}
