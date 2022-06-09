#' @title Find premature terminal codons (PTCs)
#'
#' @description Finds PTC(s) introduced by alternative exons into protein-coding transcripts.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.PSI} and \code{ParseGTF} functions.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param pval Numeric value. Adjusted p-value below which the splicing event will be analysed for PTCs.
#' @param delta Numeric value. Positive delta percent spliced-in (PSI) value above which the splicing event will be analysed for PTCs. "Positive" because only an increased in PSI value leads to increased alternative exon inclusion in the transcript.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$NMD$Prediction}.
#'
#' @importFrom Biostrings DNAStringSet translate reverseComplement
#'
#' @export


FindPTC <- function(MarvelObject, method, pval, delta) {

    # Define arguments
    method <- method
    pval <- pval
    delta <- delta

    # Example arguments
    #MarvelObject <- marvel
    #method <- c("ad", "dts")
    #pval <- c(0.10, 0.10)
    #delta <- 0.05
    
    # Subset sig. events
    .list <- list()
    
    for(i in 1:length(method)) {
        
        # Retrieve DE results table
        de.psi <- MarvelObject$DE$PSI$Table[[method[i]]]
        
        # Subset sig events
        index <- which(de.psi$p.val.adj < pval[i] & de.psi$mean.diff > delta & de.psi$outliers==FALSE)
        de.psi <- de.psi[index, ]
        
        # Subset relevant columns
        cols <- c("tran_id", "event_type", "gene_id", "gene_short_name", "gene_type")
        de.psi <- de.psi[, cols]
        
        # Save into list
        .list[[i]] <- de.psi
                
    }
    
    de.psi.small <- do.call(rbind.data.frame, .list)
    de.psi.small <- unique(de.psi.small)
        
    # Keep only SE, RI, A5SS, A3SS
    event_types <- c("SE", "RI", "A5SS", "A3SS")
    de.psi.small <- de.psi.small[which(de.psi.small$event_type %in% event_types), ]
    
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
