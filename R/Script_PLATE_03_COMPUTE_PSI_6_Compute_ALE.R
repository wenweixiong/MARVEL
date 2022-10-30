#' @title Compute alternative last exon (ALE) percent spliced-in (PSI) values
#'
#' @description Computes percent spliced-in (PSI) for alternative last exon (ALE) splicing events.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{DetectEvents} function.
#' @param CoverageThreshold Numeric value. Coverage threshold below which the PSI of the splicing event will be censored, i.e. annotated as missing (NA). Coverage defined as the total number of reads supporting both included and excluded isoforms.
#'
#' @return An object of class S3 containing with new slots \code{$SpliceFeatureValidated$ALE} and \code{$PSI$ALE}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- ComputePSI.ALE(MarvelObject=marvel.demo,
#'                               CoverageThreshold=10
#'                               )

ComputePSI.ALE <- function(MarvelObject, CoverageThreshold=10) {

    # Define arguments
    df.sj <- MarvelObject$SpliceJunction
    df.feature.posneg <- MarvelObject$SpliceFeature$ALE
    CoverageThreshold <- CoverageThreshold
    
    # Example arguments
    #MarvelObject <- marvel
    #df.sj <- MarvelObject$SpliceJunction
    #df.feature.posneg <- MarvelObject$SpliceFeature$ALE
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
        . <- strsplit(., split=":", fixed=TRUE)
        
        chr <- sapply(., function(x) {x[1]})
        start <- as.numeric(sapply(., function(x) {x[3]})) + 1
        
        # End coord
        . <- strsplit(df.feature$tran_id, split=":+@", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split="|", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split=":", fixed=TRUE)
        end <- as.numeric(sapply(., function(x) {x[2]})) - 1
        
        # Final coord
        coord.included <- paste(chr, start, end, sep=":")
        
        # Retrieve SJ counts
        sj.included <- df.sj[coord.included, ]
        
    # Retrieve excluded SJ
        # End coord
        . <- strsplit(df.feature$tran_id, split=":+@", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split="|", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split=":", fixed=TRUE)
        end <- as.numeric(sapply(., function(x) {x[1]})) - 1
        
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
    psi <- cbind.data.frame(., psi)
    row.names(psi) <- NULL
    
    # Save as new object
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
        . <- strsplit(., split="|", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split=":", fixed=TRUE)
        
        chr <- sapply(., function(x) {x[1]})
        start <- as.numeric(sapply(., function(x) {x[3]})) + 1
        
        # End coord
        . <- strsplit(df.feature$tran_id, split=":-@", fixed=TRUE)
        . <- sapply(., function(x) {x[1]})
        . <- strsplit(., split=":", fixed=TRUE)
        
        end <- as.numeric(sapply(., function(x) {x[2]})) - 1
        
        # Final coord
        coord.included <- paste(chr, start, end, sep=":")
        
        # Retrieve SJ counts
        sj.included <- df.sj[coord.included, ]
        
    # Retrieve excluded SJ
        # Start coord
        . <- strsplit(df.feature$tran_id, split=":-@", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split="|", fixed=TRUE)
        . <- sapply(., function(x) {x[2]})
        . <- strsplit(., split=":", fixed=TRUE)
        start <- as.numeric(sapply(., function(x) {x[2]})) + 1
        
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
    psi <- cbind.data.frame(., psi)
    row.names(psi) <- NULL
    
    # Save as new object
    psi.neg <- psi
    
    ######################################################################
    ###################### RETURN FINAL OBJECTS ##########################
    ######################################################################
    
    # Merge
    psi <- rbind.data.frame(psi.pos, psi.neg)
    
    # Track progress
    message(paste(nrow(psi), " splicing events validated and quantified", sep=""))
    
    # Save to new slots
    MarvelObject$SpliceFeatureValidated$ALE <- df.feature.posneg
    MarvelObject$PSI$ALE <- psi
    return(MarvelObject)

}
