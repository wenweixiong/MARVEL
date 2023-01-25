#' @title Compute skipped-exon (SE) percent spliced-in (PSI) values
#'
#' @description Validate SE splicing events and subsequently computes percent spliced-in (PSI) values these high-quality splicing events.
#'
#' @details This function computes the PSI for each SE splicing event. Splicing events provided in \code{SpliceFeature} data frame will first be cross-checked against the splice junctions provided in \code{SpliceJunction} data frame. Only events whose junctions are found in \code{SpliceJunction} are retained. The formula for computing PSI is the number of junction reads supporting the included isoform divided by the total number of reads supporting both included and excluded isoforms.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} function.
#' @param CoverageThreshold Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.
#' @param UnevenCoverageMultiplier Numeric value. Maximum allowable fold difference between two included junction counts.
#'
#' @return An object of class S3 with new slots \code{$SpliceFeatureValidated$SE}  \code{$PSI$SE}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- ComputePSI.SE(MarvelObject=marvel.demo,
#'                              CoverageThreshold=10,
#'                              UnevenCoverageMultiplier=10
#'                              )

ComputePSI.SE <- function(MarvelObject, CoverageThreshold, UnevenCoverageMultiplier=10) {

    # Define arguments
    df <- MarvelObject$SpliceFeature$SE
    sj <- MarvelObject$SpliceJunction
    CoverageThreshold <- CoverageThreshold
    UnevenCoverageMultiplier <- UnevenCoverageMultiplier
    
    # Example arguments
    #df <- marvel$SpliceFeature$SE
    #sj <- marvel$SpliceJunction
    #CoverageThreshold <- 10
    #UnevenCoverageMultiplier <- 10
    
    # Print progress
    message(paste(nrow(df), " splicing events found", sep=""))
    
    #########################################################################
    ############################# PREPARE INPUTS ############################
    #########################################################################
    
    row.names(sj) <- sj$coord.intron
    sj$coord.intron <- NULL

    # Recode missing values as no expression
    sj[is.na(sj)] <- 0
    
    #########################################################################
    ######################## FIND EVENTS IN SJ FILE #########################
    #########################################################################
    
    # +ve strand
        # Subset events
        df.pos <- df[grep(":+@", df$tran_id, fixed=TRUE), , drop=FALSE]

        # Retrieve coordinates
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        
        # Subset chr
        exon.1 <- sapply(., function(x) {x[1]})
        chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
        
        # Subset 5' included sj
        exon.1 <- sapply(., function(x) {x[1]})
        exon.2 <- sapply(., function(x) {x[2]})
        start <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[2]})) - 1
        coord.included.1 <- paste(chr, start, end, sep=":")
        
        # Subset 3' included sj
        exon.2 <- sapply(., function(x) {x[2]})
        exon.3 <- sapply(., function(x) {x[3]})
        start <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) {x[2]})) - 1
        coord.included.2 <- paste(chr, start, end, sep=":")
        
        # Subset skipped sj
        exon.1 <- sapply(., function(x) {x[1]})
        exon.3 <- sapply(., function(x) {x[3]})
        start <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) {x[2]})) - 1
        coord.excluded <- paste(chr, start, end, sep=":")
        
        # Check for sj record
        index.keep.coord.included.1 <- coord.included.1 %in% row.names(sj)
        index.keep.coord.included.2 <- coord.included.2 %in% row.names(sj)
        index.keep.coord.excluded <- coord.excluded %in% row.names(sj)
        index.keep <- (index.keep.coord.included.1 == TRUE) &
                      (index.keep.coord.included.2 == TRUE) &
                      (index.keep.coord.excluded == TRUE)
        table(index.keep.coord.included.1); table(index.keep.coord.included.2) ; table(index.keep.coord.excluded) ; table(index.keep)

        # Subset events
        df.pos <- df.pos[index.keep, , drop=FALSE]
        
    # -ve strand
        # Subset events
        df.neg <- df[grep(":-@", df$tran_id, fixed=TRUE), , drop=FALSE]

        # Retrieve coordinates
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        
        # Subset chr
        exon.1 <- sapply(., function(x) {x[1]})
        chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
        
        # Subset 5' included sj
        exon.3 <- sapply(., function(x) {x[3]})
        exon.2 <- sapply(., function(x) {x[2]})
        start <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[2]})) - 1
        coord.included.1 <- paste(chr, start, end, sep=":")
        
        # Subset 3' included sj
        exon.2 <- sapply(., function(x) {x[2]})
        exon.1 <- sapply(., function(x) {x[1]})
        start <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[2]})) - 1
        coord.included.2 <- paste(chr, start, end, sep=":")
        
        # Subset skipped sj
        exon.3 <- sapply(., function(x) {x[3]})
        exon.1 <- sapply(., function(x) {x[1]})
        start <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[2]})) - 1
        coord.excluded <- paste(chr, start, end, sep=":")
        
        # Check for sj record
        index.keep.coord.included.1 <- coord.included.1 %in% row.names(sj)
        index.keep.coord.included.2 <- coord.included.2 %in% row.names(sj)
        index.keep.coord.excluded <- coord.excluded %in% row.names(sj)
        index.keep <- (index.keep.coord.included.1 == TRUE) &
                      (index.keep.coord.included.2 == TRUE) &
                      (index.keep.coord.excluded == TRUE)
        table(index.keep.coord.included.1); table(index.keep.coord.included.2) ; table(index.keep.coord.excluded) ; table(index.keep)

        # Subset events
        df.neg <- df.neg[index.keep, , drop=FALSE]

    # Merge
    df <- rbind.data.frame(df.pos, df.neg)

    ######################################################################
    ############################ COMPUTE PSI #############################
    ######################################################################

    # +ve strand
        # Subset events
        df.pos <- df[grep(":+@", df$tran_id, fixed=TRUE), , drop=FALSE]

        # Retrieve coordinates
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        
        # Subset chr
        exon.1 <- sapply(., function(x) {x[1]})
        chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
        
        # Subset 5' included sj
        exon.1 <- sapply(., function(x) {x[1]})
        exon.2 <- sapply(., function(x) {x[2]})
        start <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[2]})) - 1
        coord.included.1 <- paste(chr, start, end, sep=":")
        
        # Subset 3' included sj
        exon.2 <- sapply(., function(x) {x[2]})
        exon.3 <- sapply(., function(x) {x[3]})
        start <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) {x[2]})) - 1
        coord.included.2 <- paste(chr, start, end, sep=":")
        
        # Subset skipped sj
        exon.1 <- sapply(., function(x) {x[1]})
        exon.3 <- sapply(., function(x) {x[3]})
        start <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) {x[2]})) - 1
        coord.excluded <- paste(chr, start, end, sep=":")

        # Retrieve sj counts
        sj.included.1 <- sj[coord.included.1, ]
        sj.included.2 <- sj[coord.included.2, ]
        sj.excluded <- sj[coord.excluded, ]

        # Compute PSI
        psi <- (sj.included.1 + sj.included.2) /
                  (sj.included.1 + sj.included.2 + (2*sj.excluded))
        psi[is.na(psi)] <- NA
        
        # Censor low coverage
        cov <- ( (sj.included.1 >= CoverageThreshold) & (sj.included.2 >= CoverageThreshold) ) |
               (sj.excluded >= CoverageThreshold)
        psi[!cov] <- NA
        
        # Censor uneven coverage
        cov <- ( (sj.included.1/sj.included.2) >= UnevenCoverageMultiplier ) |
               ( (sj.included.2/sj.included.1) >= UnevenCoverageMultiplier )
        psi[cov] <- NA
        
        # Annotate tran_id
        row.names(sj.included.1) <- df.pos$tran_id
        row.names(sj.included.2) <- df.pos$tran_id
        row.names(sj.excluded) <- df.pos$tran_id
        row.names(psi) <- df.pos$tran_id
        
        # Save as new object
        sj.included.1.pos <- sj.included.1
        sj.included.2.pos <- sj.included.2
        sj.excluded.pos <- sj.excluded
        psi.pos <- psi
        
    # -ve strand
        # Subset events
        df.neg <- df[grep(":-@", df$tran_id, fixed=TRUE), , drop=FALSE]

        # Retrieve coordinates
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        
        # Subset chr
        exon.1 <- sapply(., function(x) {x[1]})
        chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
        
        # Subset 5' included sj
        exon.3 <- sapply(., function(x) {x[3]})
        exon.2 <- sapply(., function(x) {x[2]})
        start <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[2]})) - 1
        coord.included.1 <- paste(chr, start, end, sep=":")
        
        # Subset 3' included sj
        exon.2 <- sapply(., function(x) {x[2]})
        exon.1 <- sapply(., function(x) {x[1]})
        start <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[2]})) - 1
        coord.included.2 <- paste(chr, start, end, sep=":")
        
        # Subset skipped sj
        exon.3 <- sapply(., function(x) {x[3]})
        exon.1 <- sapply(., function(x) {x[1]})
        start <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) {x[3]})) + 1
        end <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[2]})) - 1
        coord.excluded <- paste(chr, start, end, sep=":")

        # Retrieve sj counts
        sj.included.1 <- sj[coord.included.1, ]
        sj.included.2 <- sj[coord.included.2, ]
        sj.excluded <- sj[coord.excluded, ]

        # Compute PSI
        psi <- (sj.included.1 + sj.included.2) /
                  (sj.included.1 + sj.included.2 + (2*sj.excluded))
        psi[is.na(psi)] <- NA
        
        # Censor low coverage
        cov <- ( (sj.included.1 >= CoverageThreshold) & (sj.included.2 >= CoverageThreshold) ) |
               (sj.excluded >= CoverageThreshold)
        psi[!cov] <- NA
        
        # Censor uneven coverage
        cov <- ( (sj.included.1/sj.included.2) >= UnevenCoverageMultiplier ) |
               ( (sj.included.2/sj.included.1) >= UnevenCoverageMultiplier )
        psi[cov] <- NA
        
        # Annotate tran_id
        row.names(sj.included.1) <- df.neg$tran_id
        row.names(sj.included.2) <- df.neg$tran_id
        row.names(sj.excluded) <- df.neg$tran_id
        row.names(psi) <- df.neg$tran_id
        
        # Save as new object
        sj.included.1.neg <- sj.included.1
        sj.included.2.neg <- sj.included.2
        sj.excluded.neg <- sj.excluded
        psi.neg <- psi

    # Merge
    sj.included.1 <- rbind.data.frame(sj.included.1.pos, sj.included.1.neg)
    sj.included.2 <- rbind.data.frame(sj.included.2.pos, sj.included.2.neg)
    sj.excluded <- rbind.data.frame(sj.excluded.pos, sj.excluded.neg)
    psi <- rbind.data.frame(psi.pos, psi.neg)
    
    sj.included.1 <- sj.included.1[df$tran_id,]
    sj.included.2 <- sj.included.2[df$tran_id,]
    sj.excluded <- sj.excluded[df$tran_id,]
    psi <- psi[df$tran_id, ]
        
    # Remove row names
    . <- data.frame("tran_id"=row.names(sj.included.1), stringsAsFactors=FALSE)
    sj.included.1 <- cbind.data.frame(., sj.included.1)
    row.names(sj.included.1) <- NULL
    
    . <- data.frame("tran_id"=row.names(sj.included.2), stringsAsFactors=FALSE)
    sj.included.2 <- cbind.data.frame(., sj.included.2)
    row.names(sj.included.2) <- NULL
    
    . <- data.frame("tran_id"=row.names(sj.excluded), stringsAsFactors=FALSE)
    sj.excluded <- cbind.data.frame(., sj.excluded)
    row.names(sj.excluded) <- NULL
    
    . <- data.frame("tran_id"=row.names(psi), stringsAsFactors=FALSE)
    psi <- cbind.data.frame(., psi)
    row.names(psi) <- NULL
    
    # Check row orders
    table(sj.included.1$tran_id==df$tran_id)
    table(sj.included.2$tran_id==df$tran_id)
    table(sj.excluded$tran_id==df$tran_id)
    table(psi$tran_id==df$tran_id)
    
    # Print progress
    message(paste(nrow(psi), " splicing events validated and quantified", sep=""))
    
    ######################################################################
    ###################### RETURN FINAL OBJECTS ##########################
    ######################################################################
    
    # Indicate event type
    df$event_type <- "SE"
    col.others <- names(df)[-which(names(df) %in% c("tran_id", "event_type"))]
    df <- df[, c("tran_id", "event_type", col.others)]
    
    # Final files
    df.feature <- df
    df.psi <- psi
                                                    
    # Save to new slots
    MarvelObject$Counts$SE$sj.included.1 <- sj.included.1
    MarvelObject$Counts$SE$sj.included.2 <- sj.included.2
    MarvelObject$Counts$SE$sj.excluded <- sj.excluded
    MarvelObject$SpliceFeatureValidated$SE <- df.feature
    MarvelObject$PSI$SE <- psi
    return(MarvelObject)
        
}
