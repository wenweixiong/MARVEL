#' @title Compute retained-intron (RI) percent spliced-in (PSI) values
#'
#' @description Validate RI splicing events and subsequently computes percent spliced-in (PSI) values these high-quality splicing events.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param CoverageThreshold Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.
#' @param IntronCounts Data frame. Columns indicate sample IDs, rows indicate intron coordinates, and values indicate total intron coverage. The first column needs to be named \code{coord.intron}. These values will be combined with splice junction counts in the MARVEL object to compute PSI values.
#' @param thread Numeric value. Set number of threads.
#' @param read.length Numeric value. The length of read. This number will be specific to the sequencing mode. E.g. read length should be set to \code{150} when samples were sequenced in 150bp paired-end or single-end. This option should only be specified when users used read-counting approach for computing intron counts. The option should be left with its default value \code{1} when users tabulated the per-base count and summed them up to arrive at the intron counts.
#'
#' @return An object of class S3 with new slots \code{$SpliceFeatureValidated$RI} and \code{$PSI$RI}.
#'
#' @importFrom plyr join
#' @import methods
#' @importFrom utils txtProgressBar
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- ComputePSI.RI(MarvelObject=marvel.demo,
#'                              CoverageThreshold=10,
#'                              IntronCounts=marvel.demo$IntronCounts,
#'                              thread=1
#'                              )

ComputePSI.RI <- function(MarvelObject, CoverageThreshold, IntronCounts, thread, read.length=1) {

    # Define arguments
    df <- MarvelObject$SpliceFeature$RI
    sj <- MarvelObject$SpliceJunction
    CoverageThreshold <- CoverageThreshold
    df.intron.counts <- IntronCounts
    thread <- thread
    read.length <- read.length
    
    # Example arguments
    #library(parallel)
    #library(plyr)
    #df <- MarvelObject$SpliceFeature$RI
    #sj <- MarvelObject$SpliceJunction
    #CoverageThreshold <- 10
    #df.intron.counts <- MarvelObject$IntronCounts
    #thread <- 4
    #read.length <- 1
    
    if(thread==1) {
        
        message("Advisable to use >=2 threads to speed up PSI calculation")
        return(MarvelObject)
        
    }
    
    message(paste(nrow(df), " splicing events found", sep=""))

    #########################################################################
    ######################### SCALE [READ] COUNTS ###########################
    #########################################################################
    
    ## Only applicable when counts computed using read counts, instead of per-base
    
    df.intron.counts[,-1] <- df.intron.counts[,-1] * read.length
    
    #########################################################################
    ############################# PREPARE INPUTS ############################
    #########################################################################
    
    row.names(sj) <- sj$coord.intron
    sj$coord.intron <- NULL

    # Recode missing values as no expression
    sj[is.na(sj)] <- 0
    
    #########################################################################
    ###################### KEEP INDEPENDENT INTRONS #########################
    #########################################################################
    
    # Tabulate SE exons
        # Read file
        exon <- MarvelObject$SpliceFeature$SE
        # exon <- marvel$SpliceFeature$SE
        
        # +ve strand
        exon.pos <- exon[grep(":+@", exon$tran_id, fixed=TRUE), ]
        . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
        . <- unique(unlist(.))
        exon.pos <- .

        # +ve strand
        exon.neg <- exon[grep(":-@", exon$tran_id, fixed=TRUE), ]
        . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
        . <- unique(unlist(.))
        exon.neg <- .

        # Merge
        exon.se <- unique(c(exon.pos, exon.neg))
        
    # Tabulate MXE exons
        # Read file
        exon <- MarvelObject$SpliceFeature$MXE
        # exon <- marvel$SpliceFeature$MXE
        
        # +ve strand
        exon.pos <- exon[grep(":+@", exon$tran_id, fixed=TRUE), ]
        . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
        . <- unique(unlist(.))
        exon.pos <- .

        # +ve strand
        exon.neg <- exon[grep(":-@", exon$tran_id, fixed=TRUE), ]
        . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
        . <- unique(unlist(.))
        exon.neg <- .

        # Merge
        exon.mxe <- unique(c(exon.pos, exon.neg))
        
    # Tabulate RI exons
        # Read file
        exon <- df

        # +ve strand
        exon.pos <- exon[grep(":+@", exon$tran_id, fixed=TRUE), ]
        . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
        . <- unique(unlist(.))
        exon.pos <- .

        # -ve strand
            # Subset events
            exon.neg <- exon[grep(":-@", exon$tran_id, fixed=TRUE), ]
            
            # chr
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":")
            chr <- sapply(exon.1, function(x) {x[1]})

            # Exon 1
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":")
            
            start <- sapply(exon.1, function(x) {x[3]})
            
            end <- sapply(exon.1, function(x) {x[2]})
            
            exon.1 <- paste(chr, start, end, sep=":")

            # Exon 2
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            exon.2 <- strsplit(exon.2, split=":")
            
            start <- sapply(exon.2, function(x) {x[3]})
            
            end <- sapply(exon.2, function(x) {x[2]})
            
            exon.2 <- paste(chr, start, end, sep=":")

            # Merge
            exon.neg <- unique(exon.1, exon.2)

        # Merge
        exon.ri <- unique(c(exon.pos, exon.neg))

    # Tabulate A5SS exons
        # Read file
        exon <- MarvelObject$SpliceFeature$A5SS
        # exon <- marvel$SpliceFeature$A5SS

        # +ve strand
            # Subset events
            exon.pos <- exon[grep(":+@", exon$tran_id, fixed=TRUE), ]
            
            # chr
            . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":")
            chr <- sapply(exon.1, function(x) {x[1]})
            
            # Exon 1
            . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":")
            
            start <- sapply(exon.1, function(x) {x[2]})
            
            end <- sapply(exon.1, function(x) {x[3]})
            end <- strsplit(end, split="|", fixed=TRUE)
            end <- sapply(end, function(x) {x[1]})
            
            exon.1 <- paste(chr, start, end, sep=":")
            
            # Exon 2
            . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[1]})
            exon.2 <- strsplit(exon.2, split=":")
            
            start <- sapply(exon.2, function(x) {x[2]})
            
            end <- sapply(exon.2, function(x) {x[3]})
            end <- strsplit(end, split="|", fixed=TRUE)
            end <- sapply(end, function(x) {x[2]})
            
            exon.2 <- paste(chr, start, end, sep=":")
        
            # Exon 3
            . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
            exon.3 <- sapply(., function(x) {x[2]})
            
            # Merge
            exon.pos <- unique(c(exon.1, exon.2, exon.3))
            
        # -ve strand
            # Subset events
            exon.neg <- exon[grep(":-@", exon$tran_id, fixed=TRUE), ]

            # chr
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":")
            chr <- sapply(exon.1, function(x) {x[1]})

            # Exon 1
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":")
            
            end <- sapply(exon.1, function(x) {x[2]})
            
            start <- sapply(exon.1, function(x) {x[3]})
            start <- strsplit(start, split="|", fixed=TRUE)
            start <- sapply(start, function(x) {x[1]})
            
            exon.1 <- paste(chr, start, end, sep=":")
            
            # Exon 2
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[1]})
            exon.2 <- strsplit(exon.2, split=":")
            
            end <- sapply(exon.2, function(x) {x[2]})
            
            start <- sapply(exon.2, function(x) {x[3]})
            start <- strsplit(start, split="|", fixed=TRUE)
            start <- sapply(start, function(x) {x[2]})
            
            exon.2 <- paste(chr, start, end, sep=":")
            
            # Exon 3
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.3 <- sapply(., function(x) {x[2]})
            
            # Merge
            exon.neg <- unique(c(exon.1, exon.2, exon.3))
        
        # Merge
        exon.a5ss <- unique(c(exon.pos, exon.neg))

    # Tabulate A3SS exons
        # Read file
        exon <- MarvelObject$SpliceFeature$A3SS
        # exon <- marvel$SpliceFeature$A3SS

        # +ve strand
            # Subset events
            exon.pos <- exon[grep(":+@", exon$tran_id, fixed=TRUE), ]

            # chr
            . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":")
            chr <- sapply(exon.1, function(x) {x[1]})
            
            # Exon 1
            . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})

            # Exon 2
            . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})

            start <- strsplit(exon.2, split=":")
            start <- sapply(start, function(x) {x[2]})
            start <- strsplit(start, split="|", fixed=TRUE)
            start <- sapply(start, function(x) {x[1]})
            
            end <- strsplit(exon.2, split=":")
            end <- sapply(end, function(x) {x[3]})
            
            exon.2 <- paste(chr, start, end, sep=":")
            
            # Exon 3
            . <- strsplit(exon.pos$tran_id, split=":+@", fixed=TRUE)
            exon.3 <- sapply(., function(x) {x[2]})

            start <- strsplit(exon.3, split=":")
            start <- sapply(start, function(x) {x[2]})
            start <- strsplit(start, split="|", fixed=TRUE)
            start <- sapply(start, function(x) {x[2]})
            
            end <- strsplit(exon.3, split=":")
            end <- sapply(end, function(x) {x[3]})
            
            exon.3 <- paste(chr, start, end, sep=":")
        
            # Merge
            exon.pos <- unique(c(exon.1, exon.2, exon.3))
                
        # -ve strand
            # Subset events
            exon.neg <- exon[grep(":-@", exon$tran_id, fixed=TRUE), ]

            # chr
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})
            exon.1 <- strsplit(exon.1, split=":")
            chr <- sapply(exon.1, function(x) {x[1]})
            
            # Exon 1
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.1 <- sapply(., function(x) {x[1]})

            # Exon 2
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.2 <- sapply(., function(x) {x[2]})
            
            start <- strsplit(exon.2, split=":")
            start <- sapply(start, function(x) {x[3]})

            end <- strsplit(exon.2, split=":")
            end <- sapply(end, function(x) {x[2]})
            end <- strsplit(end, split="|", fixed=TRUE)
            end <- sapply(end, function(x) {x[1]})
            
            exon.2 <- paste(chr, start, end, sep=":")
            
            # Exon 3
            . <- strsplit(exon.neg$tran_id, split=":-@", fixed=TRUE)
            exon.3 <- sapply(., function(x) {x[2]})
            
            start <- strsplit(exon.3, split=":")
            start <- sapply(start, function(x) {x[3]})

            end <- strsplit(exon.3, split=":")
            end <- sapply(end, function(x) {x[2]})
            end <- strsplit(end, split="|", fixed=TRUE)
            end <- sapply(end, function(x) {x[2]})
            
            exon.3 <- paste(chr, start, end, sep=":")
            
            # Merge
            exon.neg <- unique(c(exon.1, exon.2, exon.3))
                            
        # Merge
        exon.a3ss <- unique(c(exon.pos, exon.neg))

    # Merge
    . <- unique(c(exon.se, exon.mxe, exon.ri, exon.a5ss, exon.a3ss))
    length(exon.se) ; length(exon.mxe) ; length(exon.ri) ;
     length(exon.a5ss) ; length(exon.a3ss) ; length(.)
    . <- strsplit(., split=":")
    exons.ref <- data.frame("chr"=sapply(., function(x) {x[1]}),
                            "start"=sapply(., function(x) {x[2]}),
                            "end"=sapply(., function(x) {x[3]}),
                            stringsAsFactors=FALSE
                            )
            
    # Retrieve intron coordinates for RI
        # +ve strand
            # Subset events
            df.pos <- df[grep(":+@", df$tran_id, fixed=TRUE), ]

            # Retrieve coordinates
            . <- strsplit(df.pos$tran_id, split=":+@", fixed=TRUE)
            
            # chr
            exon.1 <- sapply(., function(x) {x[1]})
            chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
            
            # Intron
            exon.1 <- sapply(., function(x) {x[1]})
            exon.2 <- sapply(., function(x) {x[2]})
            start <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[3]})) + 1
            end <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[2]})) - 1
            
            # Annotate
            df.pos$chr <- chr
            df.pos$start <- start
            df.pos$end <- end

        # -ve strand
            # Subset events
            df.neg <- df[grep(":-@", df$tran_id, fixed=TRUE), ]

            # Retrieve coordinates
            . <- strsplit(df.neg$tran_id, split=":-@", fixed=TRUE)
            
            # Subset chr
            exon.1 <- sapply(., function(x) {x[1]})
            chr <- sapply(strsplit(exon.1, ":"), function(x) {x[1]})
                
            # Subset excluded sj
            exon.2 <- sapply(., function(x) {x[2]})
            exon.1 <- sapply(., function(x) {x[1]})
            start <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) {x[2]})) + 1
            end <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) {x[3]})) - 1
            
            # Annotate
            df.neg$chr <- chr
            df.neg$start <- start
            df.neg$end <- end
            
        # Merge
        df <- rbind.data.frame(df.pos, df.neg)

    # Remove introns with overlapping exons
    #message("Retrieving independent introns... (~2mins for ~10,000 RI events)")
    
    index.keep <- NULL

    #pb <- txtProgressBar(1, nrow(df), style=3)

    for(i in 1:nrow(df)) {

        # Retrieve queries
        chr <- df[i, "chr"]
        start <- df[i, "start"]
        end <- df[i, "end"]

        # Cross-check
        exons.ref.small <- exons.ref[which(exons.ref$chr==chr), ]
        overlap <- (exons.ref.small$start >= start) &
                   (exons.ref.small$end <= end)
        overlap <- overlap[which(overlap==TRUE)]
        index.keep[i] <- ifelse(length(overlap)==0, TRUE, FALSE)
        
        # Track progress
        #setTxtProgressBar(pb, i)

    }

    table(index.keep)

    df <- df[index.keep, ]

    # Format intron coordinates
    df$coord.intron <- paste(df$chr, df$start, df$end, sep=":")
    
    #########################################################################
    ############### KEEP INDEPENDENT INTRONS (2nd PASS) #####################
    #########################################################################
    
    # Same start, different end
        # Subset events with same start
        df$chr.start <- paste(df$chr, df$start, sep=":")
        tbl <- as.data.frame(table(df$chr.start))
        coords <- as.character(tbl[which(tbl$Freq >1), "Var1"])
        
        # Split into unique and dup
        df.unique <- df[-which(df$chr.start %in% coords), ]
        df.dup <- df[which(df$chr.start %in% coords), ]
        
        if(nrow(df.dup) != 0) {
            
            # Keep longest intron for duplicates
            coords <- unique(df.dup$chr.start)
            
            .list <- list()
            
            for(i in 1:length(coords)) {
                
                # Subset duplicates
                . <- df.dup[which(df.dup$chr.start == coords[i]), ]
                
                # Keep longest
                .list[[i]] <- .[which.max(.$end),]
                
            }
            
            df.dup.collapsed <- do.call(rbind.data.frame, .list)
            
            # Merge
            df <- rbind.data.frame(df.unique, df.dup.collapsed)
            
        }

     # Same end, different start
        # Subset events with same start
        df$chr.end <- paste(df$chr, df$end, sep=":")
        tbl <- as.data.frame(table(df$chr.end))
        coords <- as.character(tbl[which(tbl$Freq >1), "Var1"])
        
        # Split into unique and dup
        df.unique <- df[-which(df$chr.end %in% coords), ]
        df.dup <- df[which(df$chr.end %in% coords), ]
        
        if(nrow(df.dup) != 0) {

            # Keep longest intron for duplicates
            coords <- unique(df.dup$chr.end)
            
            .list <- list()
            
            for(i in 1:length(coords)) {
                
                # Subset duplicates
                . <- df.dup[which(df.dup$chr.end == coords[i]), ]
                
                # Keep longest
                .list[[i]] <- .[which.min(.$start),]
                
            }
            
            df.dup.collapsed <- do.call(rbind.data.frame, .list)
            
            # Merge
            df <- rbind.data.frame(df.unique, df.dup.collapsed)
            
        }
    
    ######################################################################
    ############################ COMPUTE PSI #############################
    ######################################################################
    
    # Excluded counts
        # Create sj coord ref file
        . <- strsplit(row.names(sj), split=":")
        sj.coord <- data.frame(sapply(., function(x) {x[1]}),
                               as.numeric(sapply(., function(x) {x[2]})),
                               as.numeric(sapply(., function(x) {x[3]})),
                               stringsAsFactors=FALSE
        )
        
        names(sj.coord) <- c("chr", "start", "end")
        
        # Retrieve start, end positions
        coords <- strsplit(df$coord.intron, split=":", fixed=TRUE)
        
        # Define function
        retrieve_counts <- function(x) {
                       
            # Subset relevant junctions
            sj.coord.small <- sj.coord[which(sj.coord$chr==x[1]), ]
            overlap <- (sj.coord.small$start <= as.numeric(x[2])) &
                       (sj.coord.small$end >= as.numeric(x[3]))
            sj.coord.overlap <- sj.coord.small[overlap, ]
            sj.coord.overlap <- paste(sj.coord.overlap$chr, sj.coord.overlap$start, sj.coord.overlap$end, sep=":")
            sj.small <- sj[sj.coord.overlap, ]
            
            # Sum up counts
            counts.total <- as.data.frame(t(as.data.frame(colSums(sj.small))))
            
        }

        # Print progress
        #message("Initializing cluster with 4 nodes...")

        # Initialize cluster
        cl <- parallel::makeCluster(thread)
        parallel::clusterExport(cl=cl, varlist=c('sj.coord', 'sj'), envir=environment())
                
        # Print progress
        #message("Computing excluded counts... (~10mins for ~10,000 RI events)")
        
        # Compute average counts
        counts.total.list <- parallel::parLapply(cl, coords, retrieve_counts)

        counts.total <- do.call(rbind.data.frame, counts.total.list)
        row.names(counts.total) <- df$tran_id

        counts.excluded  <- counts.total
        
        # Stop cluster
        parallel::stopCluster(cl)
        
    # Included counts
        # Print progress
        #message("Computing included counts... (~2mins for ~10,000 RI events)")
        
        # Normalize counts by intron length
        . <- strsplit(df.intron.counts$coord.intron, split=":")
        start <- as.numeric(sapply(., function(x) {x[2]}))
        end <- as.numeric(sapply(., function(x) {x[3]}))
        length.intron <- (end - start) + 1
        
        row.names(df.intron.counts) <- df.intron.counts$coord.intron
        df.intron.counts$coord.intron <- NULL
        
        df.intron.counts.normalized <- sweep(x=df.intron.counts, MARGIN=1, FUN='/', STATS=length.intron)
        
        df.intron.counts.normalized$coord.intron <- row.names(df.intron.counts.normalized)
       
        # Annotate events
        . <- df[,c("tran_id", "coord.intron")]
        . <- join(., df.intron.counts.normalized, by="coord.intron", type="left")
        .$coord.intron <- NULL
        
        row.names(.) <- .$tran_id
        .$tran_id <- NULL
        
        # Save as new object
        counts.included <- .
        
    # Check alignment
    temp.1 <- counts.included
    temp.2 <- counts.excluded
    counts.included <- counts.included[(df$tran_id), ]
    counts.excluded <- counts.excluded[(df$tran_id), ]
    table(row.names(counts.included)==row.names(counts.excluded))
    table(row.names(counts.included)==df$tran_id)
    table(row.names(counts.excluded)==df$tran_id)
    table(names(counts.included)==names(counts.excluded))

    # Compute PSI
    #message("Computing PSI...")
    
    psi <- counts.included / (counts.included + counts.excluded)
    row.names(psi) <- df$tran_id

    # Censor psi below threshold coverage (total counts)
    cov <- counts.included + counts.excluded
    threshold <- CoverageThreshold
    psi[cov < CoverageThreshold] <- NA

    # Remove row names
    . <- data.frame("tran_id"=row.names(counts.included), stringsAsFactors=FALSE)
    counts.included <- cbind.data.frame(., counts.included)
    row.names(counts.included) <- NULL
    
    . <- data.frame("tran_id"=row.names(counts.excluded), stringsAsFactors=FALSE)
    counts.excluded <- cbind.data.frame(., counts.excluded)
    row.names(counts.excluded) <- NULL
    
    . <- data.frame("tran_id"=row.names(psi), stringsAsFactors=FALSE)
    psi <- cbind.data.frame(., psi)
    row.names(psi) <- NULL
    
    # Check row orders
    table(counts.included$tran_id==df$tran_id)
    table(counts.excluded$tran_id==df$tran_id)
    table(psi$tran_id==df$tran_id)
    
    # Print progress
    message(paste(nrow(psi), " splicing events validated and quantified", sep=""))
    
    ######################################################################
    ###################### RETURN FINAL OBJECTS ##########################
    ######################################################################
    
    # Remove intermediate columns
    df$chr <- NULL
    df$start <- NULL
    df$end <- NULL
    df$coord.intron <- NULL
    df$chr.start <- NULL
    df$chr.end <- NULL
    
    # Indicate event type
    df$event_type <- "RI"
    col.others <- names(df)[-which(names(df) %in% c("tran_id", "event_type"))]
    df <- df[, c("tran_id", "event_type", col.others)]
    
    # Final files
    df.feature <- df
    df.psi <- psi
                                        
    # Save to new slots
    MarvelObject$Counts$RI$counts.included <- counts.included
    MarvelObject$Counts$RI$counts.excluded <- counts.excluded
    MarvelObject$SpliceFeatureValidated$RI <- df.feature
    MarvelObject$PSI$RI <- psi
    
    return(MarvelObject)
    
}


