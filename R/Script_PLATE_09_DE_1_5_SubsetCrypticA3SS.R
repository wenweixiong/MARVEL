#' @title Differential gene expression analysis for differentially spliced genes
#'
#' @description Performs differential gene expression analysis between 2 groups of cells only on differentially spliced genes.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param method Vector of character string(s). To include splicing events from these method(s) for differential splicing analysis.
#' @param distance.to.ss Character string. Range of distances between A3SS and canonical splice site to consider A3SS to be cryptic. Default value \code{c(1, 100)}.
#'
#' @return An object of class S3 updated slot \code{MarvelObject$DE$PSI$Table} and new slot \code{MarvelObject$DE$PSI$A3SS.dist.to.ss}.
#'
#' @importFrom plyr join
#' @import methods
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- SubsetCrypticA3SS(MarvelObject=marvel.demo,
#'                                  method="ad",
#'                                  distance.to.ss=c(1,100)
#'                                  )
#'
#' # Check output
#' head(marvel.demo$DE$PSI$Table[["ad"]])

SubsetCrypticA3SS <- function(MarvelObject, method, distance.to.ss=c(1,100)) {

    # Define arguments
    MarvelObject <- MarvelObject
    method <- method
    distance.to.ss <- distance.to.ss
    
    # Example arguments
    #MarvelObject <- marvel
    #method <- c("ad", "dts")
    #distance.to.ss <- c(1, 100)

    ############################################################
    
    # Subset cryptic A3SS
    .list <- list()
    
    for(i in 1:length(method)) {
    
        # Subset relevent splicing DE results
        df <- MarvelObject$DE$PSI$Table[[method[i]]]
        
        # Add unique row ids for tracking
        df$row.id <- c(1:nrow(df))

        # Subset A3SS
        df.small <- df[which(df$event_type=="A3SS"), ]
        
        # Compute A3SS distance to TSS
            # +ve strand
            tran_ids.small <- grep(":+@", df.small$tran_id, fixed=TRUE, value=TRUE)
            . <- strsplit(tran_ids.small, split=":+@", fixed=TRUE)
            . <- sapply(., function(x){grep("|", x, fixed=TRUE, value=TRUE)})
            . <- strsplit(., split=":", fixed=TRUE)
            . <- sapply(., function(x){grep("|", x, fixed=TRUE, value=TRUE)})
            . <- strsplit(., split="|", fixed=TRUE)
            dist <- sapply(., function(x) {

                        abs(as.numeric(x[1]) - as.numeric(x[2]))
                    
                    })

            results.pos <- data.frame("tran_id"=tran_ids.small,
                                      "dist.to.ss"=dist,
                                      stringsAsFactors=FALSE
                                      )

            # Compute distance to canonical ss: -ve strand
            tran_ids.small <- grep(":-@", df.small$tran_id, fixed=TRUE, value=TRUE)
            . <- strsplit(tran_ids.small, split=":-@", fixed=TRUE)
            . <- sapply(., function(x){grep("|", x, fixed=TRUE, value=TRUE)})
            . <- strsplit(., split=":", fixed=TRUE)
            . <- sapply(., function(x){grep("|", x, fixed=TRUE, value=TRUE)})
            . <- strsplit(., split="|", fixed=TRUE)
            dist <- sapply(., function(x) {

                        abs(as.numeric(x[1]) - as.numeric(x[2]))
                    
                    })

            results.neg <- data.frame("tran_id"=tran_ids.small,
                                      "dist.to.ss"=dist,
                                      stringsAsFactors=FALSE
                                      )

            # Merge
            results <- rbind.data.frame(results.pos, results.neg)
            
            # Annotate splicing metadata
            results <- join(results, df.small[,c("tran_id", "event_type", "gene_id", "gene_short_name", "gene_type")], by="tran_id", type="left")
            
            cols.1 <- "dist.to.ss"
            cols.2 <- setdiff(names(results), cols.1)
            results <- results[,c(cols.2, cols.1)]
        
        # Annotate original table
        df <- join(df, results[,c("tran_id", "dist.to.ss")], by="tran_id", type="left")
        
        # Subset cryptic A3SS
        index.1 <- which(df$event_type=="A3SS" & df$dist.to.ss >= distance.to.ss[1] & df$dist.to.ss <= distance.to.ss[2])
        index.2 <- which(df$event_type != "A3SS")
        index <- c(index.1, index.2)
        df <- df[index, ]
        
        # Reorder rows
        df <- df[order(df$row.id), ]
        df$row.id <- NULL
        df$dist.to.ss <- NULL
        
        # Update slot
        MarvelObject$DE$PSI$Table[[method[i]]] <- df
        MarvelObject$DE$PSI$A3SS.dist.to.ss[[method[i]]] <- results
        
    }
    
    ############################################################
    
    # Return final object
    return(MarvelObject)
        
}
