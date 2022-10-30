#' @title Compute Alternative 3' Splice Site (A3SS) Percent Spliced-in (PSI) Values
#'
#' @description Validate A3SS splicing events and subsequently computes percent spliced-in (PSI) values these high-quality splicing events.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param CoverageThreshold Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.
#'
#' @export
#'
#' @return An object of class S3 containing with new slots \code{$SpliceFeatureValidated$A3SS} and \code{$PSI$A3SS}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- ComputePSI.A3SS(MarvelObject=marvel.demo,
#'                                CoverageThreshold=10
#'                                )

ComputePSI.A3SS <- function(MarvelObject, CoverageThreshold) {

    # Define arguments
    df <- MarvelObject$SpliceFeature$A3SS
    sj <- MarvelObject$SpliceJunction
    CoverageThreshold <- CoverageThreshold
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$SpliceFeature$A3SS
    #sj <- MarvelObject$SpliceJunction
    #CoverageThreshold <- 10
    
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

        if(nrow(df.pos) !=0) {
            
            # Subset chr
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
            
            # Subset included sj
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
            start <- as.numeric(sapply(exon.1, function(x) {x[3]})) + 1
            
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split="|", fixed=TRUE)
            exon.2 <- sapply(exon.2, function(x) {x[1]})
            exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
            end <- as.numeric(sapply(exon.2, function(x) {x[2]})) - 1

            coord.included <- paste(chr, start, end, sep=":")
                
            # Subset excluded sj
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
            start <- as.numeric(sapply(exon.1, function(x) {x[3]})) + 1
            
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split="|", fixed=TRUE)
            exon.2 <- sapply(exon.2, function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
            end <- as.numeric(sapply(exon.2, function(x) {x[1]})) - 1

            coord.excluded <- paste(chr, start, end, sep=":")
            
            # Check for sj record
            index.keep.coord.included <- coord.included %in% row.names(sj)
            index.keep.coord.excluded <- coord.excluded %in% row.names(sj)
            index.keep <- (index.keep.coord.included == TRUE) &
                          (index.keep.coord.excluded == TRUE)
            table(index.keep.coord.included); table(index.keep.coord.excluded) ; table(index.keep)

            # Subset events
            df.pos <- df.pos[index.keep, , drop=FALSE]
            
        }
        
    # -ve strand
        # Subset events
        df.neg <- df[grep(":-@", df$tran_id, fixed=TRUE), , drop=FALSE]
        
        if(nrow(df.neg) !=0) {
            
            # Subset chr
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
                
            # Subset included sj
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
            end <- as.numeric(sapply(exon.1, function(x) {x[2]})) - 1
            
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split="|", fixed=TRUE)
            exon.2 <- sapply(exon.2, function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
            start <- as.numeric(sapply(exon.2, function(x) {x[1]})) + 1

            coord.included <- paste(chr, start, end, sep=":")
                
            # Subset excluded sj
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
            end <- as.numeric(sapply(exon.1, function(x) {x[2]})) - 1
          
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split="|", fixed=TRUE)
            exon.2 <- sapply(exon.2, function(x) {x[1]})
            exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
            start <- as.numeric(sapply(exon.2, function(x) {x[2]})) + 1

            coord.excluded <- paste(chr, start, end, sep=":")
            
            # Check for sj record
            index.keep.coord.included <- coord.included %in% row.names(sj)
            index.keep.coord.excluded <- coord.excluded %in% row.names(sj)
            index.keep <- (index.keep.coord.included == TRUE) &
                          (index.keep.coord.excluded == TRUE)
            table(index.keep.coord.included); table(index.keep.coord.excluded) ; table(index.keep)

            # Subset events
            df.neg <- df.neg[index.keep, , drop=FALSE]
            
        }

    # Merge
    if(nrow(df.pos) !=0 & nrow(df.neg) !=0) {
        
        df <- rbind.data.frame(df.pos, df.neg)

    } else if(nrow(df.pos) !=0 & nrow(df.neg)==0) {
        
        df <- df.pos
        
    } else if(nrow(df.pos) == 0 & nrow(df.neg) != 0) {
        
        df <- df.neg
        
    }

    ######################################################################
    ############################ COMPUTE PSI #############################
    ######################################################################

    # +ve strand
        # Subset events
        df.pos <- df[grep(":+@", df$tran_id, fixed=TRUE), , drop=FALSE]

        if(nrow(df.pos) != 0){
            
            # Subset chr
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
            
            # Subset included sj
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
            start <- as.numeric(sapply(exon.1, function(x) {x[3]})) + 1
            
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split="|", fixed=TRUE)
            exon.2 <- sapply(exon.2, function(x) {x[1]})
            exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
            end <- as.numeric(sapply(exon.2, function(x) {x[2]})) - 1

            coord.included <- paste(chr, start, end, sep=":")
                
            # Subset excluded sj
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
            start <- as.numeric(sapply(exon.1, function(x) {x[3]})) + 1
            
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split="|", fixed=TRUE)
            exon.2 <- sapply(exon.2, function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
            end <- as.numeric(sapply(exon.2, function(x) {x[1]})) - 1

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
            
        }
        
    # -ve strand
        # Subset events
        df.neg <- df[grep(":-@", df$tran_id, fixed=TRUE), , drop=FALSE]
        
        if(nrow(df.neg) != 0){
            
            # Subset chr
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
                
            # Subset included sj
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
            end <- as.numeric(sapply(exon.1, function(x) {x[2]})) - 1
            
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split="|", fixed=TRUE)
            exon.2 <- sapply(exon.2, function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
            start <- as.numeric(sapply(exon.2, function(x) {x[1]})) + 1

            coord.included <- paste(chr, start, end, sep=":")
                
            # Subset excluded sj
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":", fixed=TRUE)
            end <- as.numeric(sapply(exon.1, function(x) {x[2]})) - 1
          
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split="|", fixed=TRUE)
            exon.2 <- sapply(exon.2, function(x) {x[1]})
            exon.2 <- strsplit(exon.2, split=":", fixed=TRUE)
            start <- as.numeric(sapply(exon.2, function(x) {x[2]})) + 1

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

        }
        
    # Merge
    if(nrow(df.pos) !=0 & nrow(df.neg) !=0) {
        
        psi <- rbind.data.frame(psi.pos, psi.neg)

    } else if(nrow(df.pos) !=0 & nrow(df.neg)==0) {
        
        psi <- psi.pos
        
    } else if(nrow(df.pos) == 0 & nrow(df.neg) != 0) {
        
        psi <- psi.neg
        
    }
    psi <- psi[df$tran_id, ]
    table(row.names(psi)==df$tran_id)

    # Remove row names
    . <- data.frame("tran_id"=row.names(psi), stringsAsFactors=FALSE)
    psi <- cbind.data.frame(., psi)
    row.names(psi) <- NULL
    
    # Print progress
    message(paste(nrow(psi), " splicing events validated and quantified", sep=""))

    ######################################################################
    ###################### RETURN FINAL OBJECTS ##########################
    ######################################################################
    
    # Indicate event type
    df$event_type <- "A3SS"
    col.others <- names(df)[-which(names(df) %in% c("tran_id", "event_type"))]
    df <- df[, c("tran_id", "event_type", col.others)]
    
    # Final files
    df.feature <- df
    df.psi <- psi
    
    # Save to new slots
    MarvelObject$SpliceFeatureValidated$A3SS <- df.feature
    MarvelObject$PSI$A3SS <- psi
    
    return(MarvelObject)
    
}


