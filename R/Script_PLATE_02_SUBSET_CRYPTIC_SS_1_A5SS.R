#' @title Subset cryptic alternative 5' splice site (A5SS)
#'
#' @description Subset A5SS splicing events whose distance to the canonical are within a user-defined distance.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CreateMarvelObject} function.
#' @param DistanceToCanonical Numeric value. Distance between alternative and canonical splice site, below which, the alternative splice site will be defined as cryptic and retained.
#'
#' @return An object of class S3 with updated slot \code{$SpliceFeature$A5SS} and new slot \code{$DistanceToCanonical$A5SS}.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export

SubsetCrypticSS.A5SS <- function(MarvelObject, DistanceToCanonical=100) {

    # Define arguments
    df.feature <- MarvelObject$SpliceFeature$A5SS
    DistanceToCanonical <- DistanceToCanonical
    
    # Example arguments
    #MarvelObject <- marvel
    #df.feature <- MarvelObject$SpliceFeature$A5SS
    #DistanceToCanonical <- 100
    
    # Add row ids for tracking
    df.feature$row_id <- c(1:nrow(df.feature))
    
    #########################################################################
    ############################ +VE STRAND #################################
    #########################################################################
    
    # Subset relevant strand
    index <- grep(":+@", df.feature$tran_id, fixed=TRUE)
    df.feature.small <- df.feature[index, ]
   
    # Retrieve ss positions
    . <- strsplit(df.feature.small$tran_id, split=":+@", fixed=TRUE)
    . <- sapply(., function(x) {x[1]})
    . <- strsplit(., split=":", fixed=TRUE)
    . <- sapply(., function(x) {x[3]})
    . <- strsplit(., split="|", fixed=TRUE)
    coord.canonical <- as.numeric(sapply(., function(x) {x[1]}))
    coord.alt <- as.numeric(sapply(., function(x) {x[2]}))
    
    # Compute distance
    dist <- abs(coord.canonical-coord.alt) + 1
    df.feature.small$dist_to_canonical <- dist
    
    # Save as new object
    df.feature.small.pos <- df.feature.small
    
    # Subset cryptic splice site
    index <- which(dist <= DistanceToCanonical)
    df.feature.small.pos. <- df.feature.small[index, ]
    
    #########################################################################
    ############################ -VE STRAND #################################
    #########################################################################
    
    # Subset relevant strand
    index <- grep(":-@", df.feature$tran_id, fixed=TRUE)
    df.feature.small <- df.feature[index, ]
   
    # Retrieve ss positions
    . <- strsplit(df.feature.small$tran_id, split=":-@", fixed=TRUE)
    . <- sapply(., function(x) {x[1]})
    . <- strsplit(., split=":", fixed=TRUE)
    . <- sapply(., function(x) {x[3]})
    . <- strsplit(., split="|", fixed=TRUE)
    coord.canonical <- as.numeric(sapply(., function(x) {x[2]}))
    coord.alt <- as.numeric(sapply(., function(x) {x[1]}))
    
    # Compute distance
    dist <- abs(coord.canonical-coord.alt) + 1
    df.feature.small$dist_to_canonical <- dist
    
    # Save as new object
    df.feature.small.neg <- df.feature.small
    
    # Subset cryptic splice site
    index <- which(dist <= DistanceToCanonical)
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
    message(paste(nrow(df.feature.small), " A5SS identified", sep=""))
    message(paste(nrow(df.feature.small.), " cryptic A5SS retained", sep=""))
    
    # Save to new slots
    MarvelObject$SpliceFeature$A5SS <- df.feature.small.
    MarvelObject$DistanceToCanonical$A5SS <- df.feature.small
    
    return(MarvelObject)

}


