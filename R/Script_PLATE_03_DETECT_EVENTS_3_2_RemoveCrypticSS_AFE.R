#' @title Remove cryptic alternative 5' splice site (A5SS)
#'
#' @description Remove cryptic A5SS splicing events from alternative first exon (AFE) splicing evnts whose distance to the canonical are within a user-defined distance.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param DistanceToCanonical Numeric value. Distance between alternative and canonical splice site, below which, the alternative splice site will be defined as cryptic and filtered out.
#'
#' @return An object of class S3 with updated slot \code{$SpliceFeature$AFE} and new slot \code{$DistanceToCanonical$AFE}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export

RemoveCrypticSS.AFE <- function(MarvelObject, DistanceToCanonical=100) {

    # Define arguments
    df.feature <- MarvelObject$SpliceFeature$AFE
    DistanceToCanonical <- DistanceToCanonical
    
    # Example arguments
    #MarvelObject <- marvel
    #df.feature <- MarvelObject$SpliceFeature$AFE
    #DistanceToCanonical <- 100
    
    # Add row ids for tracking
    df.feature$row_id <- c(1:nrow(df.feature))
    
    #########################################################################
    ############################ +VE STRAND #################################
    #########################################################################
    
    # Subset relevant strand
    index <- grep(":+@", df.feature$tran_id, fixed=TRUE)
    df.feature.small <- df.feature[index, ]
   
    # Retrieve canonical ss position
    . <- strsplit(df.feature.small$tran_id, split=":+@", fixed=TRUE)
    . <- sapply(., function(x) {x[1]})
    . <- strsplit(., split="|", fixed=TRUE)
    . <- sapply(., function(x) {x[1]})
    . <- strsplit(., split=":", fixed=TRUE)
    coord.canonical <- as.numeric(sapply(., function(x) {x[3]}))
    
    # Retrieve alt ss position
    . <- strsplit(df.feature.small$tran_id, split=":+@", fixed=TRUE)
    . <- sapply(., function(x) {x[1]})
    . <- strsplit(., split="|", fixed=TRUE)
    . <- sapply(., function(x) {x[2]})
    . <- strsplit(., split=":", fixed=TRUE)
    coord.alt <- as.numeric(sapply(., function(x) {x[2]}))
    
    # Compute distance
    dist <- abs(coord.canonical-coord.alt) + 1
    df.feature.small$dist_to_canonical <- dist
    
    # Save as new object
    df.feature.small.pos <- df.feature.small
    
    # Subset cryptic splice site
    index <- which(dist > DistanceToCanonical)
    df.feature.small.pos. <- df.feature.small[index, ]
    
    #########################################################################
    ############################ -VE STRAND #################################
    #########################################################################
    
    # Subset relevant strand
    index <- grep(":-@", df.feature$tran_id, fixed=TRUE)
    df.feature.small <- df.feature[index, ]
   
    # Retrieve canonical ss position
    . <- strsplit(df.feature.small$tran_id, split=":-@", fixed=TRUE)
    . <- sapply(., function(x) {x[1]})
    . <- strsplit(., split="|", fixed=TRUE)
    . <- sapply(., function(x) {x[2]})
    . <- strsplit(., split=":", fixed=TRUE)
    coord.canonical <- as.numeric(sapply(., function(x) {x[1]}))
    
    # Retrieve canonical alt position
    . <- strsplit(df.feature.small$tran_id, split=":-@", fixed=TRUE)
    . <- sapply(., function(x) {x[1]})
    . <- strsplit(., split="|", fixed=TRUE)
    . <- sapply(., function(x) {x[1]})
    . <- strsplit(., split=":", fixed=TRUE)
    coord.alt <- as.numeric(sapply(., function(x) {x[2]}))
    
    # Compute distance
    dist <- abs(coord.canonical-coord.alt) + 1
    df.feature.small$dist_to_canonical <- dist
    
    # Save as new object
    df.feature.small.neg <- df.feature.small
    
    # Subset cryptic splice site
    index <- which(dist > DistanceToCanonical)
    df.feature.small.neg. <- df.feature.small[index, ]
    
    #########################################################################
    ############################### MERGE ###################################
    #########################################################################
    
    # Merge filtered table
    df.feature.small. <- rbind.data.frame(df.feature.small.pos., df.feature.small.neg.)
    df.feature.small. <- df.feature.small.[order(df.feature.small.$row_id), ]
    df.feature.small.$row_id <- NULL
    df.feature.small.$dist_to_canonical <- NULL
    
    # Merge original table
    df.feature.small <- rbind.data.frame(df.feature.small.pos, df.feature.small.neg)
    df.feature.small <- df.feature.small[order(df.feature.small$row_id), ]
    df.feature.small$row_id <- NULL
    
    # Report progress
    message(paste(nrow(df.feature.small), " AFE identified", sep=""))
    message(paste(nrow(df.feature.small)-nrow(df.feature.small.), " cryptic A5SS removed", sep=""))
    
    # Save to new slots
    MarvelObject$SpliceFeature$AFE <- df.feature.small.
    MarvelObject$DistanceToCanonical$AFE <- df.feature.small
    
    return(MarvelObject)

}


