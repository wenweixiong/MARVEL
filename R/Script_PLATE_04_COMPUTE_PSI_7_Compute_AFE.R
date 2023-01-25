#' @title Compute alternative first exon (AFE) percent spliced-in (PSI) values
#'
#' @description Computes percent spliced-in (PSI) for alternative first exon (ALE) splicing events.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{DetectEvents} function.
#' @param CoverageThreshold Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.
#'
#' @return An object of class S3 containing with new slots \code{$SpliceFeatureValidated$AFE} and \code{$PSI$AFE}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- ComputePSI.AFE(MarvelObject=marvel.demo,
#'                               CoverageThreshold=10
#'                               )

ComputePSI.AFE <- function(MarvelObject, CoverageThreshold=10) {

    # Define arguments
    df.sj <- MarvelObject$SpliceJunction
    df.feature.posneg <- MarvelObject$SpliceFeature$AFE
    CoverageThreshold <- CoverageThreshold
    
    # Example arguments
    #MarvelObject <- marvel
    #df.sj <- MarvelObject$SpliceJunction
    #df.feature.posneg <- MarvelObject$SpliceFeature$AFE
    #CoverageThreshold <- 10
    
    # Print progress
    message(paste(nrow(df.feature.posneg), " splicing events found", sep=""))
    
    # Create row names
        # SJ matrix
        row.names(df.sj) <- df.sj$coord.intron
        df.sj$coord.intron <- NULL
        
    # Recode SJ NA's as 0'
    df.sj[is.na(df.sj)] <- 0

    ##########################################################################
    ############################# +VE STRAND #################################
    ##########################################################################

    # Subset relevant strand
    df.feature <- df.feature.posneg[grep(":+@", df.feature.posneg$tran_id, fixed=TRUE),]
    
    # Retrieve included SJ
        # Start coord
        . <- strsplit(df.feature$tran_id, split=":+@", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split="|", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split=":", fixed=TRUE)
        
        start <- as.numeric(sapply(., function(x) {x[2]})) + 1
        
        # End coord
        . <- strsplit(df.feature$tran_id, split=":+@", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split=":", fixed=TRUE)
        
        chr <- sapply(., function(x) {x[1]})
        end <- as.numeric(sapply(., function(x) {x[2]})) - 1
        
        # Final coord
        coord.included <- paste(chr, start, end, sep=":")
        
        # Retrieve SJ counts
        sj.included <- df.sj[coord.included, ]
        
    # Retrieve excluded SJ
        # Start coord
        . <- strsplit(df.feature$tran_id, split=":+@", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split="|", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split=":", fixed=TRUE)
        
        start <- as.numeric(sapply(., function(x) {x[3]})) + 1
        
        # Final coord
        coord.excluded <- paste(chr, start, end, sep=":")
        
        # Retrieve SJ counts
        sj.excluded <- df.sj[coord.excluded, ]
        
    # Total coverage
    sj.total <- sj.included + sj.excluded

    # Compute PSI
    psi <- sj.included/sj.total
    psi[is.na(psi)] <- NA

    # Censor low coverage
    psi[sj.total < CoverageThreshold] <- NA

    # Create tran_id column
    . <- data.frame("tran_id"=df.feature$tran_id, stringsAsFactors=FALSE)
    sj.included <- cbind.data.frame(., sj.included)
    row.names(sj.included) <- NULL
    
    . <- data.frame("tran_id"=df.feature$tran_id, stringsAsFactors=FALSE)
    sj.excluded <- cbind.data.frame(., sj.excluded)
    row.names(sj.excluded) <- NULL
    
    . <- data.frame("tran_id"=df.feature$tran_id, stringsAsFactors=FALSE)
    psi <- cbind.data.frame(., psi)
    row.names(psi) <- NULL
    
    # Save as new object
    sj.included.pos <- sj.included
    sj.excluded.pos <- sj.excluded
    psi.pos <- psi
    
    ##########################################################################
    ############################# -VE STRAND #################################
    ##########################################################################
    
    # Subset relevant strand
    df.feature <- df.feature.posneg[grep(":-@", df.feature.posneg$tran_id, fixed=TRUE),]

    # Retrieve included SJ
        # Start coord
        . <- strsplit(df.feature$tran_id, split=":-@", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split=":", fixed=TRUE)
        
        chr <- sapply(., function(x) {x[1]})
        start <- as.numeric(sapply(., function(x) {x[3]})) + 1
        
        # End coord
        . <- strsplit(df.feature$tran_id, split=":-@", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split="|", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split=":", fixed=TRUE)
        end <- as.numeric(sapply(., function(x) {x[1]})) - 1
        
        # Final coord
        coord.included <- paste(chr, start, end, sep=":")
        
        # Retrieve SJ counts
        sj.included <- df.sj[coord.included, ]
        
    # Retrieve excluded SJ
        # End coord
        . <- strsplit(df.feature$tran_id, split=":-@", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split="|", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split=":", fixed=TRUE)
        end <- as.numeric(sapply(., function(x) {x[2]})) - 1
        
        # Final coord
        coord.excluded <- paste(chr, start, end, sep=":")
        
        # Retrieve SJ counts
        sj.excluded <- df.sj[coord.excluded, ]
        
    # Total coverage
    sj.total <- sj.included + sj.excluded

    # Compute PSI
    psi <- sj.included/sj.total
    psi[is.na(psi)] <- NA

    # Censor low coverage
    psi[sj.total < CoverageThreshold] <- NA
    
    # Create tran_id column
    . <- data.frame("tran_id"=df.feature$tran_id, stringsAsFactors=FALSE)
    sj.included <- cbind.data.frame(., sj.included)
    row.names(sj.included) <- NULL
    
    . <- data.frame("tran_id"=df.feature$tran_id, stringsAsFactors=FALSE)
    sj.excluded <- cbind.data.frame(., sj.excluded)
    row.names(sj.excluded) <- NULL
    
    . <- data.frame("tran_id"=df.feature$tran_id, stringsAsFactors=FALSE)
    psi <- cbind.data.frame(., psi)
    row.names(psi) <- NULL
    
    # Save as new object
    sj.included.neg <- sj.included
    sj.excluded.neg <- sj.excluded
    psi.neg <- psi
    
    ######################################################################
    ###################### RETURN FINAL OBJECTS ##########################
    ######################################################################
    
    # Merge
    sj.included <- rbind.data.frame(sj.included.pos, sj.included.neg)
    sj.excluded <- rbind.data.frame(sj.excluded.pos, sj.excluded.neg)
    psi <- rbind.data.frame(psi.pos, psi.neg)
    
    # Check row orders
    table(sj.included$tran_id==df.feature.posneg$tran_id)
    table(sj.excluded$tran_id==df.feature.posneg$tran_id)
    table(psi$tran_id==df.feature.posneg$tran_id)
    
    # Track progress
    message(paste(nrow(psi), " splicing events validated and quantified", sep=""))
    
    # Save to new slots
    MarvelObject$Counts$AFE$sj.included <- sj.included
    MarvelObject$Counts$AFE$sj.excluded <- sj.excluded
    MarvelObject$SpliceFeatureValidated$AFE <- df.feature.posneg
    MarvelObject$PSI$AFE <- psi
    return(MarvelObject)

}
