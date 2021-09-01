#' @title Compute Alternative 5' Splice Site (A5SS) Percent Spliced-in (PSI) Values
#'
#' @description
#' \code{ComputePSI.A5SS} computes percent spliced-in (PSI) alternative 5' splice site (A5SS) splicing event.
#'
#' @details
#' This function computes the PSI for each A5SS splicing event. Splicing events provided in \code{SpliceFeature} data frame will first be cross-checked against the splice junctions provided in \code{SpliceJunction} data frame. Only events whose junctions are found in \code{SpliceJunction} are retained. The formula for computing PSI is the number of junction reads supporting the included isoform divided by the total number of reads supporting both included and excluded isoforms.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} function.
#' @param CoverageThreshold Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to two new slots. \code{$SpliceFeatureValidated$A5SS} contains the validated splicing event metadata. \code{$PSI$A5SS} contains the computed PSI values for the validated splicing events.
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @import methods
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- ComputePSI.A5SS(MarvelObject=marvel,
#'                           CoverageThreshold=10
#'                           )
#'
#' # Check output
#' marvel$SpliceFeatureValidated$A5SS
#' marvel$PSI$A5SS


ComputePSI.A5SS <- function(MarvelObject, CoverageThreshold) {

    # Define arguments
    df <- MarvelObject$SpliceFeature$A5SS
    sj <- MarvelObject$SpliceJunction
    CoverageThreshold <- CoverageThreshold
    
    print(paste("Analysing ", nrow(df), " splicing events", sep=""))
    
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

        # Subset chr
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
        
        # Subset included sj
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split="|", fixed=TRUE)
        start <- as.numeric(sapply(exon.1, function(x) {x[2]})) + 1
        
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.2 <- sapply(., function(x) {x[2]})
        exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
        end <- as.numeric(sapply(exon.2, function(x) {x[2]})) - 1

        coord.included <- paste(chr, start, end, sep=":")
            
        # Subset excluded sj
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split="|", fixed=TRUE)
        exon.1 <- sapply(exon.1, function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
        start <- as.numeric(sapply(exon.1, function(x) {x[3]})) + 1
        
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.2 <- sapply(., function(x) {x[2]})
        exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
        end <- as.numeric(sapply(exon.2, function(x) {x[2]})) - 1

        coord.excluded <- paste(chr, start, end, sep=":")
        
        # Check for sj record
        index.keep.coord.included <- coord.included %in% row.names(sj)
        index.keep.coord.excluded <- coord.excluded %in% row.names(sj)
        index.keep <- (index.keep.coord.included == TRUE) &
                      (index.keep.coord.excluded == TRUE)
        table(index.keep.coord.included); table(index.keep.coord.excluded) ; table(index.keep)

        # Subset events
        df.pos <- df.pos[index.keep, , drop=FALSE]
        
    # -ve strand
        # Subset events
        df.neg <- df[grep(":-@", df$tran_id, fixed=TRUE), , drop=FALSE]

        # Retrieve coordinates
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        
        # Subset chr
        . <- strsplit(df.neg$tran_id, split=":+@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
        
        # Subset included sj
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split="|", fixed=TRUE)
        exon.1 <- sapply(exon.1, function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
        end <- as.numeric(sapply(exon.1, function(x) {x[3]})) - 1
        
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        exon.2 <- sapply(., function(x) {x[2]})
        exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
        start <- as.numeric(sapply(exon.2, function(x) {x[3]})) + 1

        coord.included <- paste(chr, start, end, sep=":")
        
        # Subset excluded sj
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split="|", fixed=TRUE)
        end <- as.numeric(sapply(exon.1, function(x) {x[2]})) - 1

        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        exon.2 <- sapply(., function(x) {x[2]})
        exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
        start <- as.numeric(sapply(exon.2, function(x) {x[3]})) + 1

        coord.excluded <- paste(chr, start, end, sep=":")
        
        # Check for sj record
        index.keep.coord.included <- coord.included %in% row.names(sj)
        index.keep.coord.excluded <- coord.excluded %in% row.names(sj)
        index.keep <- (index.keep.coord.included == TRUE) &
                      (index.keep.coord.excluded == TRUE)
        table(index.keep.coord.included); table(index.keep.coord.excluded) ; table(index.keep)

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
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
        
        # Subset included sj
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split="|", fixed=TRUE)
        start <- as.numeric(sapply(exon.1, function(x) {x[2]})) + 1
        
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.2 <- sapply(., function(x) {x[2]})
        exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
        end <- as.numeric(sapply(exon.2, function(x) {x[2]})) - 1

        coord.included <- paste(chr, start, end, sep=":")
            
        # Subset excluded sj
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split="|", fixed=TRUE)
        exon.1 <- sapply(exon.1, function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
        start <- as.numeric(sapply(exon.1, function(x) {x[3]})) + 1
        
        . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
        exon.2 <- sapply(., function(x) {x[2]})
        exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
        end <- as.numeric(sapply(exon.2, function(x) {x[2]})) - 1

        coord.excluded <- paste(chr, start, end, sep=":")

        # Retrieve sj counts
        sj.included <- sj[coord.included, ]
        sj.excluded <- sj[coord.excluded, ]
        
        # Compute PSI
        psi <- sj.included / (sj.included + sj.excluded)
        psi[is.na(psi)] <- NA
        
        # Censor low coverage
        cov <- (sj.included >= CoverageThreshold) | (sj.excluded >= CoverageThreshold)
        psi[!cov] <- NA
        
        # Annotate tran_id
        row.names(psi) <- df.pos$tran_id
        
        # Save as new object
        psi.pos <- psi
        
    # -ve strand
        # Subset events
        df.neg <- df[grep(":-@", df$tran_id, fixed=TRUE), , drop=FALSE]

        # Retrieve coordinates
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        
        # Subset chr
        . <- strsplit(df.neg$tran_id, split=":+@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
        
        # Subset included sj
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split="|", fixed=TRUE)
        exon.1 <- sapply(exon.1, function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
        end <- as.numeric(sapply(exon.1, function(x) {x[3]})) - 1
        
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        exon.2 <- sapply(., function(x) {x[2]})
        exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
        start <- as.numeric(sapply(exon.2, function(x) {x[3]})) + 1

        coord.included <- paste(chr, start, end, sep=":")
        
        # Subset excluded sj
        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        exon.1 <- sapply(., function(x) {x[1]})
        exon.1 <- strsplit(exon.1, split="|", fixed=TRUE)
        end <- as.numeric(sapply(exon.1, function(x) {x[2]})) - 1

        . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
        exon.2 <- sapply(., function(x) {x[2]})
        exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
        start <- as.numeric(sapply(exon.2, function(x) {x[3]})) + 1

        coord.excluded <- paste(chr, start, end, sep=":")

        # Retrieve sj counts
        sj.included <- sj[coord.included, ]
        sj.excluded <- sj[coord.excluded, ]
        
        # Compute PSI
        psi <- sj.included / (sj.included + sj.excluded)
        psi[is.na(psi)] <- NA
        
        # Censor low coverage
        cov <- (sj.included >= CoverageThreshold) | (sj.excluded >= CoverageThreshold)
        psi[!cov] <- NA
        
        # Annotate tran_id
        row.names(psi) <- df.neg$tran_id
        
        # Save as new object
        psi.neg <- psi

    # Merge
    psi <- rbind.data.frame(psi.pos, psi.neg)
    psi <- psi[df$tran_id, ]
    table(row.names(psi)==df$tran_id)

    # Remove row names
    . <- data.frame("tran_id"=row.names(psi), stringsAsFactors=FALSE)
    psi <- cbind.data.frame(., psi)
    row.names(psi) <- NULL
    
    # Print progress
    print(paste(nrow(psi), " splicing events validated and quantified", sep=""))
    
    ######################################################################
    ###################### RETURN FINAL OBJECTS ##########################
    ######################################################################
    
    # Indicate event type
    df$event_type <- "A5SS"
    col.others <- names(df)[-which(names(df) %in% c("tran_id", "event_type"))]
    df <- df[, c("tran_id", "event_type", col.others)]
    
    # Final files
    df.feature <- df
    df.psi <- psi
    
    # Save to new slots
    MarvelObject$SpliceFeatureValidated$A5SS <- df.feature
    MarvelObject$PSI$A5SS <- psi
    
    return(MarvelObject)

}


